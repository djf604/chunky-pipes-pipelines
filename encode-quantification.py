import os
import subprocess
import re
import json
from datetime import datetime
from chunkypipes.components import Software, Parameter, Redirect, BasePipeline

FIRST_READS_PAIR = 0


class Pipeline(BasePipeline):
    @staticmethod
    def count_gzipped_lines(filepath):
        zcat = subprocess.Popen(['zcat', filepath], stdout=subprocess.PIPE)
        num_lines = subprocess.check_output(['wc', '-l'], stdin=zcat.stdout)
        return num_lines.strip()

    def description(self):
        return """This is an exact replication of the ENCODE long-rna pipeline.\n\n
        Requirements:\nSTAR 2.4.2a\nRSEM 1.2.15"""

    def add_pipeline_args(self, parser):
        parser.add_argument('--bam', required=True,
                            help='Path to transcriptome BAM.')
        parser.add_argument('--output', required=True,
                            help='Full path to output directory.')
        parser.add_argument('--lib', default=datetime.now().strftime('%Y-%m-%d-%H-%M-%S'),
                            help=('Name of the library, prepended to output file names. Defaults to ' +
                                  'a date string (YYYY-MM-DD-hh-mm-ss).'))
        parser.add_argument('--is-stranded', action='store_true',
                            help='Provide this argument if library is stranded.')
        parser.add_argument('--is-paired-end', action='store_true',
                            help='Provide this argument if library is paired-end.')
        return parser

    def configure(self):
        return {
            'cutadapt': {
                'path': 'Full path to cutadapt',
                'quality-base': 'Phred quality scale [33|64]'
            },
            'STAR': {
                'path': 'Full path to STAR',
                'genome-dir': 'Directory containing a STAR genome index',
                'threads': 'Number of threads to run STAR'
            },
            'RSEM': {
                'path-calculate-expression': 'Full path to RSEM (rsem-calculate-expression)',
                'path-plot-model': 'Full path to RSEM (rsem-plot-model)',
                'reference-dir': 'Full path to RSEM reference, including prefix (Ex. /path/to/rsem-ref/ref-prefix)',
                'threads': 'Number of threads to run RSEM',
                'memory': 'Memory to use for RSEM in MB'
            },
            'bedgraph_to_bw': {
                'path': 'Full path to BedgraphToBW'
            },
            'samtools': {
                'path': 'Full path to samtools'
            },
            'sort': {
                'memory': 'Memory to use for sorting (Ex. 32G)'
            },
            'bedSort': {
                'path': 'Full path to bedSort'
            }
        }

    def run_pipeline(self, pipeline_args, pipeline_config):
        # Instantiate options
        bam = pipeline_args['bam']
        output_dir = pipeline_args['output']
        logs_dir = os.path.join(output_dir, 'logs')

        # Create output, tmp, and logs directories
        subprocess.call(['mkdir', '-p', output_dir,
                         logs_dir, os.path.join(output_dir, 'tmp')])

        # Timing functions for getting running time
        start_time = datetime.now()

        # Gather QC data
        qc_data = {
            'total_raw_reads_counts': [],
            'trimmed_reads_counts': [],
            'num_reads_mapped': '0',
            'running_time_seconds': '',
            'running_time_readable': ''
        }

        # Keep list of items to delete
        staging_delete = [os.path.join(output_dir, 'tmp')]

        # Establish software instances
        rsem_calculat_expression = Software('RSEM', pipeline_config['RSEM']['path-calculate-expression'])
        rsem_plot_model = Software('RSEM', pipeline_config['RSEM']['path-plot-model'])

        # Set up RSEM parameters
        rsem_common = [
            Parameter('--bam'),
            Parameter('--estimate-rspd'),
            Parameter('--calc-ci'),
            Parameter('--no-bam-output'),
            Parameter('--seed', '12345')
        ]

        rsem_run = [
            Parameter('-p', pipeline_config['RSEM']['threads']),
            Parameter('--ci-memory', pipeline_config['RSEM']['memory'])
        ]

        rsem_type = []
        if pipeline_args['is_paired_end']:
            rsem_type.append(Parameter('--paired-end'))
        if pipeline_args['is_stranded']:
            rsem_type.append(Parameter('--forward-prob', '0'))

        # Run RSEM quantification step
        rsem_calculat_expression.run(*(rsem_common + rsem_run + rsem_type + [
            Parameter(bam),
            Parameter(pipeline_config['RSEM']['reference-dir']),
            Parameter(os.path.join(output_dir, 'RSEM_Quant')),
            Redirect(Redirect.BOTH, dest=os.path.join(logs_dir, 'Log.rsem'))
        ]))

        # Generate RSEM plot model
        rsem_plot_model.run(
            Parameter(os.path.join(output_dir, 'RSEM_Quant'), os.path.join(output_dir, 'Quant.pdf'))
        )

        # QC: Get time delta
        elapsed_time = datetime.now() - start_time
        qc_data['running_time_seconds'] = str(elapsed_time.seconds)
        qc_data['running_time_readable'] = str(elapsed_time)

        # QC: Output QC data to file
        with open(os.path.join(logs_dir, 'qc_metrics.txt'), 'w') as qc_data_file:
            qc_data_file.write(json.dumps(qc_data, indent=4) + '\n')

        # Delete temporary files
        for delete_file in staging_delete:
            subprocess.call(['rm', '-rf', delete_file])
