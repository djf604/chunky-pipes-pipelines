from chunkypipes.components import *


class Pipeline(BasePipeline):
    def description(self):
        return 'FO_NBCP filtering pipeline'

    def add_pipeline_args(self, parser):
        parser.add_argument('--mutect', help='Mutect output VCF')
        parser.add_argument('--strelka', help='Strelka output VCF')
        parser.add_argument('--dust', help='Path to dust masker file')
        parser.add_argument('--normals', action='append', help='Panel(s) of normals')
        parser.add_argument('--outfile', help='Final filtered VCF')

    def run_pipeline(self, pipeline_args, pipeline_config):
        super(Pipeline, self).run_pipeline(pipeline_args, pipeline_config)