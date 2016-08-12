from chunkypipes.components import *
# Step 0: FastQC
# Step 1a: Align with BWA, both SE and PE
# Step 1b: Post-alignment filtering

# Step 2a: Convert BAM to tagAlign, BED 3+3, both SE and PE
# Step 2b: Calculate cross-correlation QC scores
# Step 2c: Generate self-pseudoreplicates for each replicate, both SE and PE
# Step 2d: Generate pooled dataset and pooled-pseudoreplicates

# Step 3: Call peaks with MACS2, MUSIC (SPP, GEM, PeakSeq)?

# Step 4a: IDR to compare true replicates
# Step 4b: IDR to compare self-pseudoreplicates
# Step 4c: IDR to compare pooled pseudoreplicates
# Step 4d: Select final peaks calls (conservative set, optimal set)
# Step 4f: Compute IDR QC scores

# Step 6: MACS2 for histone marks
FIRST_READS_PAIR = 0


class Pipeline(BasePipeline):
    def configure(self):
        return {
            'bwa': {
                'path': 'Full path to bwa',
                'threads': 'Threads to run bwa',
                'index': 'Full path to bwa index'
            }
        }

    def description(self):
        return super(Pipeline, self).description()

    def add_pipeline_args(self, parser):
        parser.add_argument('--reads', required=True, nargs='*',
                            help='Input reads as fastq. PE separated by \':\'')
        parser.add_argument('--output', required=True,
                            help='Full path to output directory. Will be created if it doesn\'t exist.')

    def run_pipeline(self, pipeline_args, pipeline_config):
        # Get user arguments
        reads = pipeline_args['reads']
        is_paired_end = all([len(read_pair.split(':')) > 1 for read_pair in reads])

        # TODO Create output directory

        # Establish software
        bwa_aln = Software('bwa aln', pipeline_config['bwa']['path'] + ' aln')
        bwa_mem = Software('bwa mem', pipeline_config['bwa']['path'] + ' mem')

        # Step 1a: Align with bwa, both SE and PE
        for fastq in

        if is_paired_end:
            for read_pair in reads:
                for fastq in read_pair.split(':'):
                    bwa_aln.run(
                        Parameter('-q', '5'),
                        Parameter('-l', '32'),
                        Parameter('-k', '2'),
                        Parameter('-t', pipeline_config['bwa']['threads']),
                        Parameter(pipeline_config['bwa']['index']),
                        Parameter(fastq),
                        Redirect(stream=Redirect.STDOUT, dest='{}.sai'.format(fastq))
                    )
        else:
            for fastq in reads:
                bwa_aln.run(
                    Parameter('-q', '5'),
                    Parameter('-l', '32'),
                    Parameter('-k', '2'),
                    Parameter('-t', pipeline_config['bwa']['threads']),
                    Parameter(pipeline_config['bwa']['index']),
                    Parameter(fastq),
                    Redirect(stream=Redirect.STDOUT, dest='{}.sai'.format(fastq))
                )
