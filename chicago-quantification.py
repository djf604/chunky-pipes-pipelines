import os
import subprocess
import datetime
import re
import uuid
import json

from chunkypipes.components import Software, Parameter, Redirect, BasePipeline

FIRST_READS_PAIR = 0
FIRST_CHAR = 0
PERCENT_DUPLICATION = 7
MAPPED_READS_COUNT = 0

JAVA_DEFAULT_HEAP_SIZE = '6'


class Pipeline(BasePipeline):
    def description(self):
        return """RNAseq pipeline used at the University of Chicago."""

    def configure(self):
        return {
            'cufflinks': {
                'path': 'Full path to cufflinks',
                'transcriptome-gtf': 'Annotation GTF file for cufflinks',
                'threads': 'Number of threads to run cufflinks'
            },
            'htseq': {
                'path': 'Full path to htseq-count',
                'transcriptome-gtf': 'Annotation GTF file for HTSeq'
            }
        }

    def add_pipeline_args(self, parser):
        parser.add_argument('--bam', required=True,
                            help='Processed alignment BAM.')
        parser.add_argument('--output', required=True,
                            help='Full path to output directory.')
        parser.add_argument('--cufflinks-lib-type', default='fr-firststrand',
                            choices=['ff-firststrand',
                                     'ff-secondstrand',
                                     'ff-unstranded',
                                     'fr-firststrand',
                                     'fr-secondstrand',
                                     'fr-unstranded',
                                     'transfrags'],
                            help='Library type for cufflinks. Defaults to fr-firststrand.')
        parser.add_argument('--htseq-stranded', default='yes',
                            choices=['yes', 'no', 'reverse'],
                            help='Strandedness for HTSeq. Defaults to yes.')
        return parser

    def count_gzipped_lines(self, filepath):
        zcat = subprocess.Popen(['zcat', filepath], stdout=subprocess.PIPE)
        num_lines = subprocess.check_output(['wc', '-l'], stdin=zcat.stdout)
        return num_lines.strip()

    def run_pipeline(self, pipeline_args, pipeline_config):
        # Instantiate options
        bam = pipeline_args['bam']
        output_dir = pipeline_args['output']
        logs_dir = os.path.join(output_dir, 'logs')

        cufflinks_lib_type = pipeline_args['cufflinks_lib_type']
        htseq_stranded = pipeline_args['htseq_stranded']

        # Create output, tmp, and logs directories
        tmp_dir = os.path.join(output_dir, 'tmp')
        subprocess.call(['mkdir', '-p', output_dir, logs_dir, tmp_dir])

        # Keep list of items to delete
        staging_delete = [os.path.join(output_dir, 'tmp')]

        # Establish Software instances
        cufflinks = Software('Cufflinks', pipeline_config['cufflinks']['path'])
        htseq = Software('HTSeq', pipeline_config['htseq']['path'])

        cufflinks_output_dir = os.path.join(output_dir, 'cufflinks')
        subprocess.call(['mkdir', '-p', cufflinks_output_dir])
        cufflinks.run(
            Parameter('--GTF', pipeline_config['cufflinks']['transcriptome-gtf']),
            Parameter('-p', pipeline_config['cufflinks']['threads']),
            Parameter('--library-type', cufflinks_lib_type),
            Parameter('--upper-quartile-norm'),
            Parameter('-o', cufflinks_output_dir),
            Parameter('--max-bundle-frags', '1000000000'),
            Parameter(bam)
        )

        htseq_output_dir = os.path.join(output_dir, 'htseq')
        subprocess.call(['mkdir', '-p', htseq_output_dir])
        for id_attr in ['gene_id', 'gene_name']:
            for feature_type in ['gene', 'transcript', 'exon']:
                htseq.run(
                    Parameter('-f', 'bam'),
                    Parameter('-r', 'name'),
                    Parameter('-s', htseq_stranded),
                    Parameter('-t', feature_type),
                    Parameter('-i', id_attr),
                    Parameter(bam),
                    Parameter(pipeline_config['htseq']['transcriptome-gtf']),
                    Redirect(stream=Redirect.STDOUT, dest=os.path.join(htseq_output_dir,
                                                                       '{}.{}.counts'.format(feature_type,
                                                                                             id_attr)))
                )

        # Delete temporary files
        for delete_file in staging_delete:
            subprocess.call(['rm', '-rf', delete_file])
