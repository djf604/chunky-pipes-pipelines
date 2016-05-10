import os
import subprocess
from datetime import datetime
from chunkypipes.components import Software, Parameter, Redirect, BasePipeline

FIRST_READS_PAIR = 0


class Pipeline(BasePipeline):
    def description(self):
        return """RNAseq quantification pipeline using pseudo-alignments."""

    def configure(self):
        return {
            'cutadapt': {
                'path': 'Full path to cutadapt',
                'quality-base': 'Phred quality scale [33|64]',
                'minimum-length': 'Minimum length for a read to pass cutadapt'
            },
            'kallisto': {
                'path': 'Full path to kallisto',
                'index-path': 'Full path to kallisto index'
            },
            'sailfish': {
                'path': 'Full path to sailfish',
                'index-path': 'Full path to sailfish index'
            }
        }

    def add_pipeline_args(self, parser):
        parser.add_argument('--reads', required=True, action='append',
                            help=('Reads to process with this pipeline. Denote paired-end reads with ' +
                                  'a colon (Ex. read1.fastq:read2.fastq). Specify multiple times to ' +
                                  'align multiple libraries (or pairs).'))
        parser.add_argument('--output', required=True,
                            help='Full path to output directory.')
        parser.add_argument('--lib', default=datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S'),
                            help=('Name of the library, prepended to output file names. Defaults to ' +
                                  'a date string (YYYY-MM-DD-hh-mm-ss).'))
        parser.add_argument('--forward-adapter', default='ZZZ',
                            help='Adapter sequence for the forward strand.')
        parser.add_argument('--reverse-adapter', default='ZZZ',
                            help='Adapter sequnce for the reverse strand.')
        parser.add_argument('--sailfish-libtype')  # TODO Find out what the choices are
        return parser

    def run_pipeline(self, pipeline_args, pipeline_config):
        reads = pipeline_args['reads']
        output_dir = pipeline_args['output']
        logs_dir = os.path.join(output_dir, 'logs')
        lib_prefix = pipeline_args['lib']
        forward_adapter = pipeline_args['forward_adapter']
        reverse_adapter = pipeline_args['reverse_adapter']
        sailfish_libtype = pipeline_args['sailfish_libtype']

        # Determine if run is paired-end from input
        run_is_paired_end = len(reads[FIRST_READS_PAIR].split(':')) > 1

        # Create output, tmp, and logs directories
        tmp_dir = os.path.join(output_dir, 'tmp')
        subprocess.call(['mkdir', '-p', output_dir, logs_dir, tmp_dir])

        # Keep list of items to delete
        staging_delete = [os.path.join(output_dir, 'tmp')]

        cutadapt = Software('cutadapt', pipeline_config['cutadapt']['path'])
        kallisto = Software('kallisto', pipeline_config['kallisto']['path'])
        sailfish = Software('sailfish', pipeline_config['sailfish']['path'])

        # Combine reads with extra sequencing depth
        if run_is_paired_end:
            # Aggregate read1s and read2s
            read1s, read2s = [], []
            for read in reads:
                read1, read2 = read.split(':')
                read1s.append(read1)
                read2s.append(read2)

            # Combine reads groups
            combined_reads = []
            for name, reads_group in [('read1', read1s), ('read2', read2s)]:
                combined_read_filename = os.path.join(output_dir, '{}.combined.{}.fastq.gz'.format(lib_prefix, name))
                combined_reads.append(combined_read_filename)
                staging_delete.append(combined_read_filename)
                with open(combined_read_filename, 'w') as combined_reads_fastq:
                    subprocess.call(['cat'] + [read for read in reads_group], stdout=combined_reads_fastq)

            # Update reads list
            reads = ':'.join(combined_reads)
        else:
            # Combine reads
            combined_read_filename = os.path.join(output_dir, '{}.combined.fastq.gz'.format(lib_prefix))
            staging_delete.append(combined_read_filename)
            with open(combined_read_filename, 'w') as combined_reads_fastq:
                subprocess.call(['cat'] + [read for read in reads], stdout=combined_reads_fastq)

            # Update reads list
            reads = combined_read_filename

        cutadapt_common = [
            Parameter('--quality-base={}'.format(pipeline_config['cutadapt']['quality-base'])),
            Parameter('--minimum-length={}'.format(pipeline_config['cutadapt']['minimum-length'])),
            Parameter('-q', '30'),
            Redirect(stream=Redirect.STDOUT, dest=os.path.join(logs_dir, 'cutadapt.summary'))
        ]

        if run_is_paired_end:
            read1, read2 = reads.split(':')
            trimmed_read1_filename = os.path.join(output_dir, lib_prefix + '_read1.trimmed.fastq.gz')
            trimmed_read2_filename = os.path.join(output_dir, lib_prefix + '_read2.trimmed.fastq.gz')

            staging_delete.append(trimmed_read1_filename)
            staging_delete.append(trimmed_read2_filename)

            cutadapt_specific = [
                Parameter('--output={}'.format(trimmed_read1_filename)),
                Parameter('--paired-output={}'.format(trimmed_read2_filename)),
                Parameter('-a', forward_adapter),
                Parameter('-A', reverse_adapter),
                Parameter(read1),
                Parameter(read2)
            ]

            # Update reads list
            reads = ':'.join([trimmed_read1_filename, trimmed_read2_filename])
        else:
            # Construct new filename
            trimmed_read_filename = os.path.join(output_dir, lib_prefix + '.trimmed.fastq.gz')

            staging_delete.append(trimmed_read_filename)

            cutadapt_specific = [
                Parameter('--output={}'.format(trimmed_read_filename)),
                Parameter('-a', forward_adapter),
                Parameter(reads)
            ]

            # Update reads list
            reads = [trimmed_read_filename]

        # Run cutadapt
        cutadapt.run(*(cutadapt_common + cutadapt_specific))

        # Step 3: Kallisto Quantification
        kallisto_common = [
            Parameter('--index={}'.format(pipeline_config['kallisto']['index-path'])),
            Parameter('--output-dir={}'.format(os.path.join(output_dir, 'kallisto_quant')))
        ]

        if run_is_paired_end:
            read1, read2 = reads.split(':')
            kallisto_ended = [
                Parameter(read1),
                Parameter(read2)
            ]
        else:
            kallisto_ended = [
                Parameter(reads)
            ]

        # Run kallisto
        kallisto.run(*(kallisto_common + kallisto_ended))

        # Step 4: Sailfish Quantification
        sailfish_common = [
            Parameter('--index', pipeline_config['sailfish']['index-path']),
            Parameter('--libType', '"{}"'.format(sailfish_libtype)),
            Parameter('--output', os.path.join(output_dir, 'sailfish_quant'))
        ]

        if run_is_paired_end:
            read1, read2 = reads.split(':')
            sailfish_ended = [
                Parameter('-1', '<(zcat {})'.format(read1)),
                Parameter('-2', '<(zcat {})'.format(read2)),
            ]
        else:
            sailfish_ended = [
                Parameter('-r', '<(zcat {})'.format(reads))
            ]

        sailfish.run(*(sailfish_common + sailfish_ended), shell=True)

        # Delete staged items
        for item in staging_delete:
            subprocess.call(['rm', '-rf', item])









# def run_pipeline(reads, options):
#     # Instantiate options
#     output_dir = options['output_dir']
#     logs_dir = options['logs_dir']
#     lib_prefix = options['lib_prefix']
#     step = options['step']
#     config = options['config']
#     run_is_paired_end = options['run_is_paired_end']
#
#     # Keep list of items to delete
#     staging_delete = ['tmp']
#
#     try:
#         forward_adapter = options['extra_info']['forward_adapter']
#         reverse_adapter = options['extra_info']['reverse_adapter']
#         sailfish_libtype = options['extra_info']['sailfish_libtype']
#     except KeyError, e:
#         # TODO Make better exception message
#         raise KeyError('Some needed extra_info not given.')
#
#     # Establish Software instances
#     cat = Software('cat', '/bin/cat')
#     cutadapt = Software('cutadapt', config['cutadapt']['path'])
#     kallisto = Software('kallisto', config['kallisto']['path'])
#     sailfish = Software('sailfish', config['sailfish']['path'])
#
#     # Combine reads with extra sequencing depth
#     if step <= 1 and len(reads) >= 2:
#         if run_is_paired_end:
#             # Aggregate read1s and read2s
#             read1s, read2s = [], []
#             for read in reads:
#                 read1, read2 = read.split(':')
#                 read1s.append(read1)
#                 read2s.append(read2)
#
#             # Combine reads groups
#             combined_reads = []
#             for name, reads_group in [('read1', read1s), ('read2', read2s)]:
#                 combined_read_filename = os.path.join(output_dir, '{}.combined.{}.fastq.gz'.format(lib_prefix, name))
#                 combined_reads.append(combined_read_filename)
#                 staging_delete.append(combined_read_filename)
#                 cat.run(
#                     Parameter(*[read for read in reads_group]),
#                     Redirect(type='1>', dest=combined_read_filename)
#                 )
#
#             # Update reads list
#             reads = [':'.join(combined_reads)]
#         else:
#             # Combine reads
#             combined_read_filename = os.path.join(output_dir, '{}.combined.fastq.gz'.format(lib_prefix))
#             staging_delete.append(combined_read_filename)
#             cat.run(
#                 Parameter(*[read for read in reads]),
#                 Redirect(type='1>', dest=combined_read_filename)
#             )
#
#             # Update reads list
#             reads = [combined_read_filename]
#
#     # Trim adapters with cutadapt
#     if step <= 2:
#         reads = reads[0]
#         if run_is_paired_end:
#             # Get paired-end reads, construct new filenames
#             read1, read2 = reads.split(':')
#             trimmed_read1_filename = os.path.join(output_dir, lib_prefix + '_read1.trimmed.fastq.gz')
#             trimmed_read2_filename = os.path.join(output_dir, lib_prefix + '_read2.trimmed.fastq.gz')
#
#             staging_delete.append(trimmed_read1_filename)
#             staging_delete.append(trimmed_read2_filename)
#
#             # Run cutadapt
#             cutadapt.run(
#                 Parameter('--quality-base={}'.format(config['cutadapt']['quality-base'])),
#                 Parameter('--minimum-length=5'),
#                 Parameter('--output={}'.format(trimmed_read1_filename)),
#                 Parameter('--paired-output={}'.format(trimmed_read2_filename)),
#                 Parameter('-a', forward_adapter if forward_adapter else 'ZZZ'),
#                 Parameter('-A', reverse_adapter if reverse_adapter else 'ZZZ'),
#                 Parameter('-q', '30'),
#                 Parameter(read1),
#                 Parameter(read2),
#                 Redirect(type='1>', dest=os.path.join(output_dir, 'logs', 'cutadapt.trendy.summary'))
#             )
#
#             # Update reads list
#             reads = ':'.join([trimmed_read1_filename, trimmed_read2_filename])
#
#         else:
#             # Construct new filename
#             trimmed_read_filename = os.path.join(output_dir, lib_prefix + '.trimmed.fastq.gz')
#
#             staging_delete.append(trimmed_read_filename)
#
#             # Run cutadapt
#             cutadapt.run(
#                 Parameter('--quality-base={}'.format(config['cutadapt']['quality-base'])),
#                 Parameter('--minimum-length=5'),
#                 Parameter('--output={}'.format(trimmed_read_filename)),
#                 Parameter('-a', forward_adapter if forward_adapter else 'ZZZ'),
#                 Parameter('-q', '30'),
#                 Parameter(reads[0]),
#                 Redirect(type='1>', dest=os.path.join(output_dir, 'logs', 'cutadapt.trendy.summary'))
#             )
#
#             # Update reads list
#             reads = [trimmed_read_filename]
#
#     # Run Kallisto
#     if step <= 3:
#         if run_is_paired_end:
#             read1, read2 = reads.split(':')
#             kallisto.run(
#                 Parameter('--index={}'.format(config['kallisto']['index-path'])),
#                 Parameter('--output-dir={}'.format('kallisto_quant')),
#                 Parameter(read1),
#                 Parameter(read2)
#             )
#         else:
#             kallisto.run(
#                 Parameter('--index={}'.format(config['kallisto']['index-path'])),
#                 Parameter('--output-dir={}'.format('kallisto_quant')),
#                 Parameter(reads[0])
#             )
#
#     # Run Sailfish
#     if step <= 4:
#         if run_is_paired_end:
#             read1, read2 = reads.split(':')
#             sailfish.run(
#                 Parameter('--index', config['sailfish']['index-path']),
#                 Parameter('--libType', '\"{}\"'.format(sailfish_libtype)),
#                 Parameter('-1', '<(zcat {})'.format(read1)),
#                 Parameter('-2', '<(zcat {})'.format(read2)),
#                 Parameter('--output', 'sailfish_quant')
#             )
#         else:
#             sailfish.run(
#                 Parameter('--index', config['sailfish']['index-path']),
#                 Parameter('--libType', '\"{}\"'.format(sailfish_libtype)),
#                 Parameter('-r', '<(zcat {})'.format(reads[0])),
#                 Parameter('--output', 'sailfish_quant')
#             )
#
#         # Delete staged items
#         for item in staging_delete:
#             subprocess.call(['rm', '-rf', item])
