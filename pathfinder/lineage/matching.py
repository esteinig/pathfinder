"""
===================================
Pathfinder Lineage Matching Module
===================================

Lineage matching prototype. Tries to determine lineage based on nanopore
reads and a database of MRSA lineages or other `pf-survey` results.
"""

from old.process.results import SurveyResult

from pathlib import Path
import pandas

from colorama import Fore

Y = Fore.YELLOW
R = Fore.RED
LR = Fore.LIGHTRED_EX
G = Fore.GREEN
C = Fore.CYAN
CY = Fore.LIGHTCYAN_EX
RE = Fore.RESET
MA = Fore.LIGHTMAGENTA_EX
YE = Fore.LIGHTYELLOW_EX
BL = Fore.LIGHTBLUE_EX
LG = Fore.LIGHTGREEN_EX


class LineageMatching:

    def __init__(
            self,
            survey: SurveyResult = None
    ):
        self.survey = survey

    def read_survey(
            self,
            survey_result: Path
    ) -> None:

        self.survey = SurveyResult(
            path=Path()
        )
        self.survey.data.read(
            survey_result
        )

    def select_sequence_types(
        self,
        outdir: Path = Path.home() / 'transfer' / 'lmdb' / 'seqs',
        min_count=None,
        sample=None,
        values=None,
    ):

        data = self.survey.data.select(
            'mlst', 'sequence_type', values=values,
            min_count=min_count, sample=sample
        )

        mlst = data.mlst.sequence_type.sort_index()
        pheno = data.mykrobe_phenotype.sort_index()

        sep = pandas.Series(
            ['_' for _ in data.iid.iid], index=mlst.index
        )

        index = mlst.index.to_series()

        index += sep + mlst.astype(str) + sep
        for column in pheno.columns:
            index += pheno[column]

        data.link_fasta(fdir=outdir, index=index.to_dict(), symlink=True)
