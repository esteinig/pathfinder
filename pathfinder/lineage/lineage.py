"""
==========================
Pathfinder Lineage Module
==========================

Lineage analyser prototype.
"""

from old.process.results import SurveyResult
from pathfinder.database.client import PathfinderClient


class LineageAnalyser(PathfinderClient):
    """ Lineage specific summaries and analyses on `pandas.DataFrame` objects
    or directly from the `pathfinder.database.client.Client` for MongoDB

    Read the data object that inherits from a pipeline ResultData as in
     `pathfinder.process.data.SurveyData`. These are simply dataclass
     containers with `pandas.DataFrame` attributes that are filled in
     the parse an processing and cleaning steps or read from storage.
    """

    def __init__(self, survey):

        self.survey: SurveyResult = survey

        self._group_data_mlst()

    def _group_data_mlst(self):

        """ Group all dataframes in the data object by sequence type """

        self.survey.data.groupby('mlst', 'sequence_type', set_index=True)

    def _summarize_genomes(self):

        pass
