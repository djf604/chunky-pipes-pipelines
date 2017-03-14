import os
import pysam
import subprocess
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
            'fastqc': {
                'path': 'Full path to FastQC'
            },
            'picard': {
                'path': 'Full path to Picard',
                'CollectRnaSeqMetrics': {
                    'ref-flat': 'Full path to refFlat file',
                    'ribosomal-intervals': 'Full path to ribosomal intervals file'
                },
            },
            'preseq': {
                'path': 'Full path to preseq',
                'bam2mr': 'Full path to bam2mr'
            },
            'RNA-SeQC': {
                'path': 'Full path to Broad RNA-SeQC'
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

        fastqc_output_dir = os.path.join(pipeline_args['output_dir'], 'fastqc')
        subprocess.call('mkdir -p {}'.format(fastqc_output_dir), shell=True)
        for fastq in pipeline_args['fastqs']:
            fastqc.run(
                Parameter('--outdir={}'.format(fastqc_output_dir)),
                Parameter('--threads', '8'),
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
            - CollectGcBiasMetrics
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
        subprocess.call('mkdir -p {}'.format(picard_output_dir), shell=True)

        # Create interval list
        tmp_interval_list = os.path.join(pipeline_config['picard']['CollectRnaSeqMetrics']['ribosomal-intervals'] +
                                         '.tmp.intervals')
        subprocess.call('samtools view -H {} > header.tmp.txt'.format(sorted_bam), shell=True)
        subprocess.call('cat header.tmp.txt {} > {}'.format(
            pipeline_config['picard']['CollectRnaSeqMetrics']['ribosomal-intervals'],
            tmp_interval_list
        ), shell=True)
        os.remove('header.tmp.txt')

        picard['MarkDuplicates'].run(
            Parameter('INPUT={}'.format(sorted_bam)),
            Parameter('OUTPUT=/dev/null'),
            Parameter('METRICS_FILE={}'.format(os.path.join(picard_output_dir, 'markduplicates.metrics'))),
            Parameter('TMP_DIR=/mnt/analysis/tmp')
        )

        picard['CollectRnaSeqMetrics'].run(
            Parameter('REF_FLAT={}'.format(pipeline_config['picard']['CollectRnaSeqMetrics']['ref-flat'])),
            Parameter('RIBOSOMAL_INTERVALS={}'.format(tmp_interval_list)),
            Parameter('STRAND_SPECIFICITY={}'.format(
                'SECOND_READ_TRANSCRIPTION_STRAND' if pipeline_args['is_stranded'] else 'NONE'
            )),
            Parameter('INPUT={}'.format(sorted_bam)),
            Parameter('OUTPUT={}'.format(os.path.join(picard_output_dir, 'rnaseq.metrics'))),
            Parameter('TMP_DIR=/mnt/analysis/tmp')
        )

        picard['CollectInsertSizeMetrics'].run(
            Parameter('INPUT={}'.format(sorted_bam)),
            Parameter('OUTPUT={}'.format(os.path.join(picard_output_dir, 'insert_size.metrics'))),
            Parameter('HISTOGRAM_FILE=/dev/null'),
            Parameter('TMP_DIR=/mnt/analysis/tmp')
        )

        picard['CollectAlignmentSummaryMetrics'].run(
            Parameter('REFERENCE_SEQUENCE={}'.format(pipeline_config['reference-genome'])),
            Parameter('INPUT={}'.format(sorted_bam)),
            Parameter('OUTPUT={}'.format(os.path.join(picard_output_dir, 'alignment_summary.metrics'))),
            Parameter('TMP_DIR=/mnt/analysis/tmp')
        )

        picard['CollectGcBiasMetrics'].run(
            Parameter('INPUT={}'.format(sorted_bam)),
            Parameter('OUTPUT={}'.format(os.path.join(picard_output_dir, 'gcbias.metrics'))),
            Parameter('CHART_OUTPUT=/dev/null'),
            Parameter('SUMMARY_OUTPUT={}'.format(os.path.join(picard_output_dir, 'gcbias.summary'))),
            Parameter('REFERENCE_SEQUENCE={}'.format(pipeline_config['reference-genome'])),
            Parameter('TMP_DIR=/mnt/analysis/tmp')
        )

        picard['EstimateLibraryComplexity'].run(
            Parameter('INPUT={}'.format(sorted_bam)),
            Parameter('OUTPUT={}'.format(os.path.join(picard_output_dir, 'library_complexity.metrics'))),
            Parameter('TMP_DIR=/mnt/analysis/tmp')
        )

        os.remove(tmp_interval_list)

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
        picard = kwargs['picard']
        samtools_faidx = kwargs['samtools_faidx']
        pipeline_config = kwargs['pipeline_config']
        pipeline_args = kwargs['pipeline_args']
        sorted_bam = kwargs['sorted_bam']

        # If reference hasn't been indexed, run samtools faidx to index
        if not os.path.isfile(pipeline_config['reference-genome'] + '.fai'):
            samtools_faidx.run(Parameter(pipeline_config['reference-genome']))

        # If reference doesn't have a dictionary files, run Picard CreateSequenceDictionary
        if not os.path.isfile(pipeline_config['reference-genome'] + '.dict'):
            picard['CreateSequenceDictionary'].run(
                Parameter('REFERENCE={}'.format(pipeline_config['reference-genome'])),
                Parameter('OUTPUT={}'.format(pipeline_config['reference-genome'] + '.dict'))
            )

        # Add read groups to BAM
        tmp_readgroups_bam = sorted_bam + '.tmp.readgroups.bam'
        picard['AddOrReplaceReadGroups'].run(
            Parameter('INPUT={}'.format(sorted_bam)),
            Parameter('OUTPUT={}'.format(tmp_readgroups_bam)),
            Parameter('RGLB={}'.format(pipeline_args['lib'])),
            Parameter('RGPL=illumina'),
            Parameter('RGPU=flowcellid'),
            Parameter('TMP_DIR=/mnt/analysis/tmp'),
            Parameter('RGSM={}'.format(pipeline_args['lib']))
        )
        subprocess.call('samtools index {}'.format(tmp_readgroups_bam), shell=True)

        # Run RNA-SeQC
        rnaseqc_output_dir = os.path.join(pipeline_args['output_dir'], 'RNA-SeQC')
        subprocess.call('mkdir -p {}'.format(rnaseqc_output_dir), shell=True)
        rnaseqc.run(
            Parameter('-o', rnaseqc_output_dir),
            Parameter('-r', pipeline_config['reference-genome']),
            Parameter('-s', '"{sample_id}|{bam_path}|None"'.format(
                sample_id=pipeline_args['lib'],
                bam_path=tmp_readgroups_bam
            )),
            Parameter('-t', pipeline_config['transcriptome-gtf']),
            Parameter('-singleEnd') if not pipeline_args['is_paired_end'] else Parameter()
        )

        os.remove(tmp_readgroups_bam)
        os.remove(tmp_readgroups_bam + '.bai')

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
        subprocess.call('mkdir -p {}'.format(preseq_output_dir), shell=True)

        preseq['c_curve'].run(
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

        featurecounts_output_dir = os.path.join(pipeline_args['output_dir'], 'featurecounts.txt')
        subprocess.call('mkdir -p {}'.format(featurecounts_output_dir), shell=True)
        featurecounts.run(
            Parameter('-a', pipeline_config['transcriptome-gtf']),  # Annotation file
            Parameter('-o', featurecounts_output_dir),  # Output file
            Parameter('-s', '2') if pipeline_args['is_stranded'] else Parameter(),
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
            in {'CreateSequenceDictionary', 'MarkDuplicates', 'CollectRnaSeqMetrics',
                'CollectInsertSizeMetrics', 'CollectAlignmentSummaryMetrics', 'CollectGcBiasMetrics',
                'EstimateLibraryComplexity', 'AddOrReplaceReadGroups'}
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

        # Create output directory
        subprocess.call('mkdir -p {}'.format(pipeline_args['output_dir']), shell=True)
        subprocess.call('mkdir -p /mnt/analysis/tmp', shell=True)

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
            picard=picard,
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

        # self.run_preseq(
        #     preseq=preseq,
        #     bam2mr=bam2mr,
        #     sorted_bam=sorted_bam,
        #     pipeline_args=pipeline_args
        # )

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
        subprocess.call('rm -rf /mnt/analysis/tmp', shell=True)