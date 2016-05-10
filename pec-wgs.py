import os
from chunkypipes.components import *


class Pipeline(BasePipeline):
    def configure(self):
        return {
            'cutadapt': {
                'path': 'Full path to cutadapt'
            },
            'bwa': {
                'path': 'Full path to bwa',
                'index-prefix': 'Full path to bwa index prefix'
            },
            'samtools': {
                'path': 'Full path to samtools'
            }
        }

    def add_pipeline_args(self, parser):
        parser.add_argument('--forward-adapter')
        parser.add_argument('--reverse-adapter')
        parser.add_argument('--reads')
        parser.add_argument('--output')
        parser.add_argument('--lib')

    def description(self):
        return 'Using BWA to align PEC Whole-Genome Sequence data'

    def run_pipeline(self, pipeline_args, pipeline_config):
        reads = pipeline_args['reads']
        forward_adapter = pipeline_args['forward_adapter']
        reverse_adapter = pipeline_args['reverse_adapter']
        output_dir = os.path.abspath(pipeline_args['output'])
        lib_prefix = pipeline_args['lib']

        logs_dir = os.path.join(output_dir, 'logs')

        subprocess.call(['mkdir', '-p', output_dir, logs_dir])

        cutadapt = Software('cutadapt', pipeline_config['cutadapt']['path'])
        bwa_mem = Software('bwa mem', pipeline_config['bwa']['path'] + ' mem')
        samtools_view = Software('samtools', pipeline_config['samtools']['path'] + 'view')

        trimmed_fastq_filenames = [os.path.join(output_dir, lib_prefix + '.trimmed.R{}.fastq.gz')
                                   for i, _ in enumerate(reads)]
        cutadapt.run(
            Parameter('--minimum-length=5'),
            Parameter('-q', '30'),
            Parameter('--output={}'.format(trimmed_fastq_filenames[0])),
            Parameter('--paired-output={}'.format(trimmed_fastq_filenames[1])),
            Parameter('-a', 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'),
            Parameter('-A', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA'),
            Parameter(reads[0]),
            Parameter(reads[1]),
            Redirect(stream=Redirect.STDOUT, dest=os.path.join(logs_dir, 'cutadapt.summary.log'))
        )

        bwa_mem.run(
            Parameter('-t', pipeline_config['bwa']['threads']),
            Parameter(pipeline_config['bwa']['index-prefix']),
            Parameter(trimmed_fastq_filenames[0]),
            Parameter(trimmed_fastq_filenames[1]),
            Redirect(stream=Redirect.STDERR, dest=os.path.join(logs_dir, 'bwa_mem.log')),
            Pipe(
                samtools_view.pipe(
                    Parameter('-hSb'),
                    Parameter('-o', 'output.bam'),
                    Parameter('-')
                )
            )
        )
