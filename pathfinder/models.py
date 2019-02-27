from mongoengine import *

import os
import tarfile

from pathlib import Path
from pathfinder.process.parsers import get_content, PipelineProcessor

import logging

logging.basicConfig(level=logging.INFO, filename='model.log', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

from tqdm import tqdm

class Pipeline:

    def __init__(self, result_path):

        self.name = "survey"
        self.version = "0.2"
        self.result_path = result_path

    def decompress(self, path:str) -> Path:

        path = Path(path)
        results = Path(str(self.result_path))
        if str(results).endswith('tar.gz'):
            with tarfile.open(results, 'r:gz') as tar:
                tar.extractall(path=path)

        return path

    def _prep_parse(self, pipeline_processor, decompress: bool,
                    tmpdir: str) -> (PipelineProcessor, dict):

        if pipeline_processor is None:
            pipeline_processor = PipelineProcessor()

        if decompress:
            os.makedirs(tmpdir, exist_ok=True)
            path = self.decompress(path=tmpdir)
        else:
            path = self.result_path

        dirs = get_content(path=path)

        return pipeline_processor, dirs


class SurveyModel(Pipeline):

    def process(self, pipeline_processor=None, decompress: bool=False, tmpdir: str='tmp'):

        processes = ('mlst', 'kraken', 'mash', 'resfinder', 'vfdb', 'plasmidfinder', 'mykrobe_susceptibility',
                     'mykrobe_genotype')

        pp, dirs = self._prep_parse(pipeline_processor, decompress, tmpdir)

        results = {}
        # Aggregate = True for all parsers
        for process, files in tqdm(dirs.items()):
            if process.name == 'mlst':
                results[process.name] = pp.parse(files=files, remove='.tab', parse_func=pp.parse_mlst)
            elif process.name == 'kraken':
                results[process.name] = pp.parse(files=files, remove='.report', parse_func=pp.parse_kraken,
                                                 process_func=pp.process_kraken)
            elif process.name == 'mash':
                results[process.name] = pp.parse(files=files, remove='.mash.tab', parse_func=pp.parse_mash)
            elif process.name == 'mykrobe':

                # TODO - multiple outputs from parse function into aggregation?
                results['mykrobe_genotype'] = pp.parse(files=files, remove='.json',
                                                       parse_func=pp.parse_mykrobe_genotype)

                results['mykrobe_susceptibility'] = pp.parse(files=files, remove='.json',
                                                             parse_func=pp.parse_mykrobe_susceptibility)

            elif process.name == 'resfinder':
                # Directly accessing subdirectory for database:
                results[process.name] = pp.parse(files, remove='.tab', parse_func=pp.parse_abricate,
                                                 process_func=pp.process_abricate)
            elif process.name == 'vfdb':
                # Directly accessing subdirectory for database:
                results[process.name] = pp.parse(files=files, remove='.tab', parse_func=pp.parse_abricate,
                                                 process_func=pp.process_abricate)
            elif process.name == 'plasmidfinder':
                # Directly accessing subdirectory for database:
                results[process.name] = pp.parse(files=files, remove='.tab', parse_func=pp.parse_abricate,
                                                 process_func=pp.process_abricate)
            else:
                pass

        # Checking if all results present:
        for process in processes:
            if process not in results.keys():
                results[process] = None
                print(f'Could not detect results from {process} in Survey.')

        return results

    def write_results(self, results, outdir):

        outpath = Path(outdir)
        outpath.mkdir(exist_ok=True, parents=True)

        for process, result in results.items():
            result.to_csv(f'{outpath / process}.csv', index=False)




    def create_genomes(self, results: dict):

        ids = self._get_unique_indices(results)

        #
        # for id in ids:
        #     genome = Genome(id=id, accession=id)
        #     genome.parse_survey(results=results)

    def create_survey(self):
        pass

    @staticmethod
    def _get_unique_indices(results: dict) -> list:

        unique_ids = []
        for process, result in results.items():
            for id in result.index.unique():
                if id not in unique_ids:
                    unique_ids.append(id)

        return unique_ids

# class Genome(Document):
#
#     _id = ObjectIdField()
#     _survey = ObjectIdField()
#
#     id = StringField()
#     accession = StringField()
#
#     mlst = StringField()
#     species = StringField()
#
#     genotype = EmbeddedDocument(Genotype)
#
#     assembly = DictField()
#     stats = DictField()
#     mash = DictField()
#
#     # Data Fields for more thorough,
#     # but unstructured data (DictFields)
#
#     qc = DictField()
#     data = DictField()
#     reads = DictField()
#
#     def parse_survey(self, results):
#
#         data = {
#             'kraken': self.from_kraken(results.get('kraken', None), level="S", top=3),
#             'mlst': self.from_mlst(results.get('mlst', None)),
#             'resfinder': None,
#             'vfdb': None,
#             'plasmidfinder': None,
#             'mykrobe': None,
#             'fastqc_before': None,
#             'fastqc_after': None,
#             'mash': None,
#             'prokka': None
#         }
#
#     def from_kraken(self, df, level="S", top=3):
#
#         if df is None or df.empty:
#             pass
#
#     def from_mlst(self, df):
#
#         if df is not None:
#             # Transform to dictionary:
#             mlst = df.drop(["file", "species"], axis=1).to_dict(orient="index")[self.id]
#             # Assign sequence type to Genome
#             self.mlst = mlst["sequence_type"]
#             return mlst
#
#         else:
#             self.mlst = None
#             return None
#
#     def from_abricate(self):
#
#         pass
#
#
#     def _check_and_subset_data(self, attr, df):
#
#         pass
#
#     def _compute_stats(self):
#         pass
#
#     def _plot_summary(self):
#         pass
#
#     def _plot_genotypes(self):
#         pass
#
#     def _plot_stats(self):
#         pass
