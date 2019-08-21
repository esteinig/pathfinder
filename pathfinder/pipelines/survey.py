from pathlib import Path
from pathfinder.pipelines.data import SurveyData
from pathfinder.utils import get_subdict, retain_files, get_id_from_fname
from tqdm import tqdm

import dataclasses
import pandas
import json

from dataclasses import dataclass
from typing import Callable


@dataclass
class SurveyProcessSetting:
    """ Template for processe pipeline data """

    files: list = None
    parse_func: Callable = None
    process_func: Callable = None
    remove: str = '.extension'


class SurveyResult:

    """ Parse results from: pf-core/pf-survey """

    def __init__(self, path: Path = None):

        self.path = path
        self.data: SurveyData = SurveyData()

    def __add__(self, other):

        # Enable sum() with __radd__
        # by skipping starting condition:
        if other == 0:
            return self
        else:
            if not isinstance(other, SurveyResult):
                raise TypeError('Only other SurveyResult objects can be added.')
            self.path = None
            self.data += other.data

            return self

    def __radd__(self, other):

        return self.__add__(other)

    # Main access methods

    def parse(
        self
    ):
        sp = SurveyParser(path=self.path)

        for process, df in sp.parse():
            setattr(self.data, process, df)

        # self.data.update_iid()

        # print(self.data)

    def filter(
        self,
        species: int or str = None,
        purity: float = 0.8,
        contamination: float = 0.02,
        gene_identity: float = 0.9,
        gene_coverage: float = 0.9,
        complete_lineage: bool = True
    ):

        sf = SurveyFilter(
            species=species,
            purity=purity,
            contamination=contamination,
            gene_identity=gene_identity,
            gene_coverage=gene_coverage,
            complete_lineage=complete_lineage
        )

        self.data.remove(
            dict(
                kraken=sf.clean_kraken(self.data.kraken),
                mlst=sf.clean_mlst(self.data.mlst)
            )
        )


#################
# Survey Filter #
#################

@dataclasses.dataclass
class SurveyFilter:

    purity: float = 0.8
    species: int or str = None
    contamination: float = 0.02
    gene_identity: float = 0.9
    gene_coverage: float = 0.9
    complete_lineage: bool = True

    def clean_kraken(
            self, data: pandas.DataFrame
    ) -> list:

        """ Clean and condense the taxonomic read assignments by Kraken2

        Inspects each isolate's results from the agglomerated data frame
        and identify the following:


            .. py:attr:`pathfinder.pipelines.SurveyFilter.contamination`

                Threshold for removing contamination with this percentage of
                reads classified in top species besides the most represented
                species.

            .. py:attr:`pathfinder.pipelines.SurveyFilter.purity`

                Threshold for retaining samples with at least this percentage
                of reads assigned to the top species.

        :returns A list of isolate identifiers to be removed at the
            end of the cleaning process.

        """

        isolates = data.groupby(by=data.index)

        uncertain, contaminated, misidentified = 0, 0, 0,

        to_remove = list()
        for _, g in isolates:

            if not g.empty:
                iid = g.index[0]

                g = g.loc[g['level'] == "S"]

                g = g.nlargest(3, 'percent')

                if g['percent'].iloc[0] < self.purity * 100:
                    to_remove.append(iid)
                    uncertain += 1

                if (g['percent'][1:] > self.contamination * 100).any():
                    contaminated += 1
                    to_remove.append(iid)

                if self.species:
                    if isinstance(self.species, str):
                        column = 'taxonomy'
                    else:
                        column = 'taxid'

                    if g[column].iloc[0] != self.species:
                        to_remove.append(iid)
                        misidentified += 1

        return list(
            pandas.Series(to_remove).unique()
        )

    def clean_mlst(self, data: pandas.DataFrame) -> list:

        if self.complete_lineage:
            return data[data['sequence_type'] == '-'].index.tolist()
        else:
            return list()


#################
# Survey Parser #
#################


class SurveyParser:
    """ Collection of parsing methods and data frame headers
    to make life easier by using an instance of the SurveyParser """

    def __init__(self, path: Path):

        # Files to parse

        self.path = path

        # DataFrame columns: Surveys

        self.mlst_columns = [
            "file", "species", "sequence_type",
        ]
        self.mash_columns = [
            "file", "ref", "dist", "p-value", "match",
        ]
        self.kraken_report_columns = [
            "percent", "reads", "direct", "level", "taxid", "taxonomy",
        ]
        self.prokka_columns = [
            "locus_tag", "feature", "length", "gene", "ec", "cog", "product",
        ]
        self.abricate_columns = [
            "file", "sequence", "start", "end", "gene", "coverage_bases",
            "coverage_map", "gaps",  "coverage", "identity", "database",
            "accession", "product",
        ]
        self.kleborate_columns = [
            "strain", "species", "st", "virulence_score", "resistance_score",
            "Yersiniabactin", "YbST", "Colibactin", "CbST", "Aerobactin",
            "AbST", "Salmochelin", "SmST", "rmpA", "rmpA2", "wzi",  "K_locus",
            "K_locus_confidence", "O_locus", "O_locus_confidence", "AGly",
            "Col", "Fcyn", "Flq", "Gly", "MLS", "Ntmdz", "Phe", "Rif", "Sul",
            "Tet", "Tmt", "Bla", "Bla_Carb", "Bla_ESBL", "Bla_ESBL_inhR",
            "Bla_broad", "Bla_broad_inhR",
        ]

    # Main access parsing method:

    def parse(self):

        """ Parse - Process - Clean """

        processes = {

            'mykrobe_lineage': SurveyProcessSetting(
                files=list(
                    (self.path / 'mykrobe').glob('*.json')
                ),
                parse_func=self.parse_mykrobe_lineage,
                remove='.json'
            ),
            'kraken': SurveyProcessSetting(
                files=list(
                    (self.path / 'kraken').glob('*.report')
                ),
                parse_func=self.parse_kraken,
                process_func=self.process_kraken,
                remove='.report'
            ),
            'mlst': SurveyProcessSetting(
                files=list(
                    (self.path / 'mlst').glob('*.tab')
                ),
                parse_func=self.parse_mlst,
                remove='.tab'
            ),
            'mash': SurveyProcessSetting(
                files=list(
                    (self.path / 'mash').glob('*.mash.tab')
                ),
                parse_func=self.parse_mash,
                remove='.mash.tab'
            ),
            'mykrobe_phenotype': SurveyProcessSetting(
                files=list(
                    (self.path / 'mykrobe').glob('*.json')
                ),
                parse_func=self.parse_mykrobe_susceptibility,
                remove='.json'
            ),
            'mykrobe_genotype': SurveyProcessSetting(
                files=list(
                    (self.path / 'mykrobe').glob('*.json')
                ),
                parse_func=self.parse_mykrobe_genotype,
                remove='.json'
            ),
            'abricate_resistance': SurveyProcessSetting(
                files=list(
                    (self.path / 'abricate' / 'resfinder').glob('*.tab')
                ),
                parse_func=self.parse_abricate,
                process_func=self.process_abricate,
                remove='.tab'
            ),
            'abricate_virulence': SurveyProcessSetting(
                files=list(
                    (self.path / 'abricate' / 'vfdb').glob('*.tab')
                ),
                parse_func=self.parse_abricate,
                process_func=self.process_abricate,
                remove='.tab'
            ),
            'abricate_plasmid': SurveyProcessSetting(
                files=list(
                    (self.path / 'abricate' / 'plasmidfinder').glob('*.tab')
                ),
                parse_func=self.parse_abricate,
                process_func=self.process_abricate,
                remove='.tab'
            ),
            'kleborate': SurveyProcessSetting(
                files=list(
                    (self.path / 'kleborate').glob('*.tab')
                ),
                parse_func=self.parse_kleborate,
                remove='.tab'
            )
        }

        for process, process_setting in processes.items():
            try:
                yield process, self._parse_and_process(
                    **dataclasses.asdict(process_setting),
                    process=process
                )
            except ValueError:
                print(f'Could not process: {process}')
                continue

    # Private helper methods for parsing

    def _parse_and_process(
        self,
        files: [],
        parse_func: Callable = None,
        process_func: Callable or None = None,
        remove: str or list = None,
        subset: str = None,
        process: str = None,
        *args, **kwargs
    ) -> dict or pandas.DataFrame:

        """Anonymous function to parse results based on a list of files,
        the name components to be stripped for extracting IDs, and the
        parse method of this class to return a dictionary of:

        :param files:
            list of input files from glob for parsing.
        :param remove:
            str or list of str, remove from filename for ID.
        :param parse_func:
            parse method of this class for parsing.
        :param process_func:
            process of data frame after reading
        :param subset:
            parse only files that contain substring

        :returns parsed results (dict) or aggregated results (dataframe)
        """

        files = list(files)

        if subset:
            files = retain_files(files, retain=subset)

        # Parse results from many files by ID:
        data = {
            get_id_from_fname(file, remove=remove): parse_func(file, *args)
            for file in tqdm(files, desc=f'{process}')
        }

        # Aggregate the data into DataFrames:
        data = self._aggregate(data)

        if process_func:
            data = process_func(data, **kwargs)

        return data

    def _aggregate(self, parsed_data):
        """ Aggregate results from the parsed data"""

        ids, dfs = self._split_data(parsed_data)

        df = pandas.concat(
            [
                df.assign(
                    id=[ids[i] for _ in range(len(df))]
                )
                for i, df in enumerate(dfs)
            ], sort=False).set_index('id')

        df.index.name = "ID"

        return df

    @staticmethod
    def _split_data(data) -> iter:
        """Split a data dictionary (keys: ids, values: data)
         into two id-sorted lists """

        return zip(*[(ide, data[ide]) for ide in sorted(data)])

    # Public helper methods for result file parsing

    def parse_kraken(self, file) -> pandas.DataFrame:
        """Parse single report output from Kraken"""
        return pandas.read_csv(
            file, header=None, sep="\t", names=self.kraken_report_columns
        )

    def parse_mash(self, file) -> pandas.DataFrame:
        """Parse single output from Mash"""
        return pandas.read_csv(
            file, sep="\t", header=None, names=self.mash_columns
        )

    def parse_abricate(self, file) -> pandas.DataFrame:
        """Parse single output from Abricate"""
        return pandas.read_csv(
            file, sep="\t", header=0, names=self.abricate_columns
        )

    def parse_prokka(self, file, gff=True) -> pandas.DataFrame:
        """Parse single output feature table from Prokka"""

        if gff:
            return self._parse_prokka_gff(file)
        else:
            return pandas.read_csv(
                file, sep="\t", header=0, names=self.prokka_columns
            )

    @staticmethod
    def _parse_prokka_gff(file):

        data = []
        with open(file, 'r') as infile:
            for line in infile:
                if line.startswith('##'):
                    pass
                elif line.startswith('##FASTA'):
                    return 0
                else:
                    data.append(
                        line.strip().split('\t')
                    )

        return pandas.DataFrame(data, columns=[
            'contig', 'inference', 'type', 'start', 'stop', 'unknown_1',
            'strand', 'unknown_2', 'description'
        ])

    def parse_mlst(self, file) -> pandas.DataFrame:
        """Parse single output from mlst (github.com/tseemann/mlst)"""

        df = pandas.read_csv(file, sep="\t", header=None)

        df.columns = self.mlst_columns + [
            str(i) for i, col in enumerate(df.columns, 1)
        ][:-3]

        return df

    def parse_kleborate(self, file) -> pandas.DataFrame:
        """ Parse single file output from Kleborate """

        return pandas.read_csv(
            file, sep='\t', header=0
        )

    def parse_mykrobe_susceptibility(self, file) -> pandas.DataFrame:

        susceptibilities, _ = self._parse_mykrobe(file)

        return susceptibilities

    def parse_mykrobe_genotype(self, file) -> pandas.DataFrame:

        _, genotype = self._parse_mykrobe(file)

        return genotype

    def parse_mykrobe_lineage(self, file) -> pandas.DataFrame:

        return self._parse_mykrobe(file, lineage=True)

    @staticmethod
    def _parse_mykrobe(
            file, lineage=False
    ) -> (pandas.DataFrame, pandas.DataFrame) or pandas.DataFrame:
        """Parse single output from MykrobePredictor, only parse
        called resistances and the genes that called them.
        """

        with open(file, "r") as infile:
            mykrobe = json.load(infile)

            if lineage:
                # TODO make sure order highest if multiple?
                keys = [
                    k for k, v in list(
                        get_subdict("lineage", mykrobe)
                    )[0].items()
                ]
                return pandas.DataFrame.from_dict({'lineage': keys})
            else:
                susceptibility = list(
                    get_subdict("susceptibility", mykrobe)
                )[0]

                data = {
                    drug: {
                        "susceptibility": data["predict"],
                        "genotype": list(
                            data["called_by"].keys()
                        )
                        if "called_by" in data.keys() else None
                    }
                    for drug, data in susceptibility.items()
                }

                susceptibilities = pandas.DataFrame.from_dict(
                    {
                        drug: list(dat["susceptibility"])
                        for drug, dat in data.items()
                    }
                )

                genotypes = pandas.DataFrame.from_dict(
                    {
                        drug: [None] if not dat["genotype"]
                        else [",".join(dat["genotype"])]
                        for drug, dat in data.items()
                     }
                )

                return susceptibilities, genotypes

    # Public methods for processing

    @staticmethod
    def process_abricate(df):
        # Remove alleles from gene names:
        df.gene = df.gene.str.rsplit('_').str[0]
        return df

    @staticmethod
    def process_kraken(df):
        df.taxonomy = df.taxonomy.str.strip()
        return df