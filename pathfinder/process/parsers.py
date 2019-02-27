import json
import pandas

from pathfinder.utils import get_subdict


class Parser:
    """ Just a collection of parsing functions and data frame headers
    to make life easier by using an instance of the Parser class """

    def __init__(self):

        # DataFrame columns: Surveys

        self.mlst_columns = [
            "file", "species", "sequence_type",
        ]
        self.mash_columns = [
            "file", "ref", "dist", "p-value", "match",
        ]
        self.fastqc_summary_columns = [
            "alert", "category", "file",
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

    #####################################
    # Parse methods for SurveyProcessor #
    #####################################

    def parse_fastqc(self, file, summary=True) -> pandas.DataFrame:
        """ Parse FastQC output, currently only combined
        (cat) summary files for forward and reverse reads """

        if summary:
            return pandas.read_csv(
                file, sep="\t", header=None, names=self.fastqc_summary_columns
            )

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

    def parse_mykrobe_susceptibility(self, file) -> pandas.DataFrame:

        susceptibilities, _ = self._parse_mykrobe(file)

        return susceptibilities

    def parse_mykrobe_genotype(self, file) -> pandas.DataFrame:

        _, genotype = self._parse_mykrobe(file)

        return genotype

    @staticmethod
    def _parse_mykrobe(file) -> (pandas.DataFrame, pandas.DataFrame):
        """Parse single output from MykrobePredictor, only parse
        called resistances and the genes that called them.
        """

        with open(file, "r") as infile:
            mykrobe = json.load(infile)
            susceptibility = list(get_subdict("susceptibility", mykrobe))[0]

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

