import os
import pysam
from chunkypipes.components import Software, Parameter, BasePipeline


class Pipeline(BasePipeline):
    def description(self):
        return 'QC pipeline for the PsychENCODE Data Analysis Core'

    def add_pipeline_args(self, parser):
        parser.add_argument('--bam', help='Full path to BAM file')
        parser.add_argument('--fastqs', nargs='*', help='Space separated list of full paths to fastq files')
        parser.add_argument('--lib', help='Sample name')
        parser.add_argument('--output-dir', help='Root output directory for QC programs')
        parser.add_argument('--is-paired-end', action='store_true', help='Whether sample was sequenced as paired-end')
        parser.add_argument('--is-stranded', action='store_true', help=('Whether library was created with a stranded '
                                                                        'protocol'))

    def configure(self):
        return {
            'fastq': {
                'path': 'Full path to FastQC'
            },
            'picard': {
                'path': 'Full path to Picard',
                'MarkDuplicates': {
                    'path': 'Full path to Picard MarkDuplicates'
                },
                'CollectRnaSeqMetrics': {
                    'path': 'Full path to Picard CollectRnaSeqMetrics',
                    'ref-flat': 'Full path to refFlat file',
                    'ribosomal-intervals': 'Full path to ribosomal intervals file'
                },
                'CollectInsertSizeMetrics': {
                    'path': 'Full path to Picard CollectInsertSizeMetrics'
                },
                'CollectAlignmentSummaryMetrics': {
                    'path': 'Full path to Picard CollectAlignmentSummaryMetrics'
                },
                'CollectGCBiasMetrics': {
                    'path': 'Full path to Picard CollectGCBiasMetrics'
                },
                'EstimateLibraryComplexity': {
                    'path': 'Full path to Picard EstimateLibraryComplexity'
                }
            },
            'preseq': {
                'path': 'Full path to preseq',
                'bam2mr': 'Full path to bam2mr'
            },
            'RNA-SeQC': {
                'path': 'Full path to Broad RNA-SeQC',
                'reference': 'Full path to reference genome in fasta format',
                'transcripts': 'Full path to GTF file defining transcripts'
            },
            'featureCounts': {
                'path': 'Full path to featureCounts'
            },
            'samtools': {
                'path': 'Full path to samtools'
            },
            'novosort': {
                'path': 'Full path to novosort'
            },
            'reference-genome': 'Full path to reference genome fasta',
            'transcriptome-gtf': 'Full path to GTF file defining transcripts'
        }

    @staticmethod
    def run_fastqc(**kwargs):
        """
        Run FastQC over all fastqs associated with this sample.
        Expects the following keyword arguments:
            - fastq::Software
            - pipeline_args
        """
        fastqc = kwargs['fastqc']
        pipeline_args = kwargs['pipeline_args']

        for fastq in pipeline_args['fastqs']:
            fastqc.run(
                Parameter('--outdir={}'.format(os.path.join(pipeline_args['output_dir'], 'fastqc'))),
                Parameter(fastq)
            )

    @staticmethod
    def run_picard_suite(**kwargs):
        """
        Run a desired subset of the Picard suite:
            - MarkDuplicates
            - CollectRnaSeqMetrics
            TODO get RefFlat file, ribosomal intervals
            - CollectInsertSizeMetrics
            - CollectAlignmentSummaryMetrics
            - CollectGCBiasMetrics
            - EstimateLibraryComplexity
        Expects the following keyword arguments:
            - picard::dict<Software> a dictionary where each key is the name of the subprogram
                                     and the value is the Software instance
            - sorted_bam::str
            - pipeline_config
            - pipeline_args
        """
        picard = kwargs['picard']
        sorted_bam = kwargs['sorted_bam']
        pipeline_config = kwargs['pipeline_config']
        pipeline_args = kwargs['pipeline_args']

        picard_output_dir = os.path.join(pipeline_args['output_dir'], 'picard')

        picard['MarkDuplicates'].run(
            Parameter('INPUT={}'.format(sorted_bam)),
            Parameter('OUTPUT=/dev/null'),
            Parameter('METRICS_FILE={}'.format(os.path.join(picard_output_dir, 'markduplicates.metrics')))
        )

        picard['CollectRnaSeqMetrics'].run(
            Parameter('REF_FLAT={}'.format(pipeline_config['picard']['CollectRnaSeqMetrics']['ref-flat'])),
            Parameter('RIBOSOMAL_INTERVALS={}'.format(
                pipeline_config['picard']['CollectRnaSeqMetrics']['ribosomal-intervals']
            )),
            Parameter('STRAND_SPECIFICITY={}'.format(
                'SECOND_READ_TRANSCRIPTION_STRAND' if pipeline_args['is_stranded'] else 'NONE'
            )),
            Parameter('INPUT={}'.format(sorted_bam)),
            Parameter('OUTPUT={}'.format(os.path.join(picard_output_dir, 'rnaseq.metrics')))
        )

        picard['CollectInsertSizeMetrics'].run(
            Parameter('INPUT={}'.format(sorted_bam)),
            Parameter('OUTPUT={}'.format(os.path.join(picard_output_dir, 'insert_size.metrics'))),
            Parameter('HISTOGRAM_FILE=/dev/null'),
        )

        picard['CollectAlignmentSummaryMetrics'].run(
            Parameter('REFERENCE_SEQUENCE={}'.format(pipeline_config['reference-genome'])),
            Parameter('INPUT={}'.format(sorted_bam)),
            Parameter('OUTPUT={}'.format(os.path.join(picard_output_dir, 'alignment_summary.metrics')))
        )

        picard['CollectGCBiasMetrics'].run(
            Parameter('INPUT={}'.format(sorted_bam)),
            Parameter('OUTPUT={}'.format(os.path.join(picard_output_dir, 'gcbias.metrics'))),
            Parameter('CHART_OUTPUT=/dev/null'),
            Parameter('SUMMARY_OUTPUT={}'.format(os.path.join(picard_output_dir, 'gcbias.summary')))
        )

        picard['EstimateLibraryComplexity'].run(
            Parameter('INPUT={}'.format(sorted_bam)),
            Parameter('OUTPUT={}'.format(os.path.join(picard_output_dir, 'library_complexity.metrics')))
        )


    @staticmethod
    def run_rnaseqc(**kwargs):
        """
        Run the Broad program RNA-SeQC, which has the worst name on the planet.
        It's a java jar, but other than that doesn't have any dependencies.
        It has a lot of pre-requisities about the BAM and reference files, which are taken care of here
        Expects the follow keyword arguments:
            - rnaseqc::Software
            - picard_csd::Software
            - samtools_index::Software
            - samtools_faidx::Software
            - pipeline_config
            - pipeline_args
            - sorted_bam::str
        """
        rnaseqc = kwargs['rnaseqc']
        picard_csd = kwargs['picard_csd']
        samtools_faidx = kwargs['samtools_faidx']
        pipeline_config = kwargs['pipeline_config']
        pipeline_args = kwargs['pipeline_args']
        sorted_bam = kwargs['sorted_bam']

        # If reference hasn't been indexed, run samtools faidx to index
        if not os.path.isfile(pipeline_config['RNA-SeQC']['reference'] + '.fai'):
            samtools_faidx.run(Parameter(pipeline_config['RNA-SeQC']['reference']))

        # If reference doesn't have a dictionary files, run Picard CreateSequenceDictionary
        if not os.path.isfile(pipeline_config['RNA-SeQC']['reference'] + '.dict'):
            picard_csd.run(
                Parameter('REFERENCE={}'.format(pipeline_config['RNA-SeQC']['reference'])),
                Parameter('OUTPUT={}'.format(pipeline_config['RNA-SeQC']['reference'] + '.dict')),
                Parameter('NUM_SEQUENCES=null')
            )

        # Run RNA-SeQC
        rnaseqc.run(
            Parameter('-o', os.path.join(pipeline_args['output_dir'], 'RNA-SeQC')),
            Parameter('-r', pipeline_config['reference-genome']),
            Parameter('-s', '"{sample_id}|{bam_path}|None"'.format(
                sample_id=pipeline_args['lib'],
                bam_path=sorted_bam
            )),
            Parameter('-t', pipeline_config['transcriptome-gtf']),
            Parameter('-singleEnd') if not pipeline_args['is_paired_end'] else Parameter()
        )

    @staticmethod
    def run_preseq(**kwargs):
        """
        Runs the preseq program
        Expects the following keyword arguments:
            - preseq::Software
            - bam2mr::Software
            - sorted_bam::str
            - pipeline_config
            - pipeline_args
        """
        preseq = kwargs['preseq']
        bam2mr = kwargs['bam2mr']
        sorted_bam = kwargs['sorted_bam']
        pipeline_args = kwargs['pipeline_args']

        preseq_output_dir = os.path.join(pipeline_args['output_dir'], 'preseq')
        preseq['c_count'].run(
            Parameter('-output', os.path.join(preseq_output_dir, 'c_count.txt')),
            Parameter('-bam'),
            Parameter('-pe') if pipeline_args['is_paired_end'] else Parameter(),
            Parameter(sorted_bam)
        )

        preseq['lc_extrap'].run(
            Parameter('-output', os.path.join(preseq_output_dir, 'lc_extrap.txt')),
            Parameter('-bam'),
            Parameter('-pe') if pipeline_args['is_paired_end'] else Parameter(),
            Parameter(sorted_bam)
        )

        mr_file = os.path.join(pipeline_args['output_dir'], 'tmp.mr')
        bam2mr.run(
            Parameter('-output', mr_file),
            Parameter(sorted_bam)
        )
        preseq['gc_extrap'].run(
            Parameter('-output', os.path.join(preseq_output_dir, 'gc_extrap.txt')),
            Parameter(mr_file)
        )
        os.remove(mr_file)

    @staticmethod
    def run_featurecounts(**kwargs):
        """
        Runs the featureCounts program, stranded and paired-end if needed
        Expects the following keyword arguments:
            - featurecounts::Software
            - pipeline_args
            - pipeline_config
        """
        featurecounts = kwargs['featurecounts']
        sorted_bam = kwargs['sorted_bam']
        pipeline_args = kwargs['pipeline_args']
        pipeline_config = kwargs['pipeline_config']

        featurecounts.run(
            Parameter('-a', pipeline_config['transcriptome-gtf']),  # Annotation file
            Parameter('-o', os.path.join(pipeline_args['output_dir'], 'featurecounts.txt')),  # Output file
            Parameter('-s 2') if pipeline_args['is_stranded'] else Parameter(),
            Parameter('-p') if pipeline_args['is_paired_end'] else Parameter(),
            Parameter(sorted_bam)
        )

    @staticmethod
    def run_chrm_percentage(**kwargs):
        sorted_bam = kwargs['sorted_bam']
        pipeline_args = kwargs['pipeline_args']
        sorted_pybam = pysam.AlignmentFile(sorted_bam, 'rb')

        with open(os.path.join(pipeline_args['output_dir'], 'chrM.txt'), 'w') as chrm:
            chrm.write('{}\n'.format(sorted_pybam.count(reference='chrM')/float(sorted_pybam.mapped)))

    def run_pipeline(self, pipeline_args, pipeline_config):
        # Instantiate Software instances
        fastqc = Software('FastQC', pipeline_config['fastqc']['path'])
        rnaseqc = Software('RNA-SeQC', pipeline_config['RNA-SeQC']['path'])

        picard = {
            subprogram_name: Software('picard {}'.format(subprogram_name),
                                      pipeline_config['picard']['path'] + ' {}'.format(subprogram_name))
            for subprogram_name
            in {'CollectSequenceDictionary', 'MarkDuplicates', 'CollectRnaSeqMetrics',
                'CollectInsertSizeMetrics', 'CollectAlignmentSummaryMetrics', 'CollectGCBiasMetrics',
                'EstimateLibraryComplexity'}
        }

        preseq = {
            subprogram_name: Software('preseq {}'.format(subprogram_name),
                                      pipeline_config['preseq']['path'] + ' {}'.format(subprogram_name))
            for subprogram_name
            in {'c_curve', 'lc_extrap', 'gc_extrap'}
        }
        bam2mr = Software('bam2mr', pipeline_config['preseq']['bam2mr'])

        featurecounts = Software('featureCounts', pipeline_config['featureCounts']['path'])

        samtools_faidx = Software('samtools faidx', pipeline_config['samtools']['path'] + ' faidx')
        novosort = Software('novosort', pipeline_config['novosort']['path'])

        # Sort bam file
        sorted_bam = os.path.join(pipeline_args['output_dir'], 'sorted.tmp.bam')
        novosort.run(
            Parameter('--index'),
            Parameter('--output', sorted_bam),
            Parameter(pipeline_args['bam'])
        )

        # Run FastQC
        self.run_fastqc(
            fastqc=fastqc,
            pipeline_args=pipeline_args
        )

        # Run RNA-SeQC
        self.run_rnaseqc(
            rnaseqc=rnaseqc,
            picard_csd=picard['CollectSequenceDictionary'],
            samtools_faidx=samtools_faidx,
            pipeline_config=pipeline_config,
            pipeline_args=pipeline_args,
            sorted_bam=sorted_bam
        )

        # Run Picard suite
        self.run_picard_suite(
            picard=picard,
            sorted_bam=sorted_bam,
            pipeline_config=pipeline_config,
            pipeline_args=pipeline_args
        )

        self.run_preseq(
            preseq=preseq,
            bam2mr=bam2mr,
            sorted_bam=sorted_bam,
            pipeline_args=pipeline_args
        )

        self.run_featurecounts(
            featurecounts=featurecounts,
            sorted_bam=sorted_bam,
            pipeline_args=pipeline_args,
            pipeline_config=pipeline_config
        )

        self.run_chrm_percentage(
            sorted_bam=sorted_bam,
            pipeline_args=pipeline_args
        )

        # Remove temporary sorted bam
        os.remove(sorted_bam)
        os.remove(sorted_bam + '.bai')