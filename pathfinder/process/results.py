"""
==========================
Pathfinder Results Module
==========================

Processes output of Nextflow pipelines into standardized format
and links to the database models.

"""

import pandas
from pathlib import Path

from pathfinder.process.processors import SurveyProcessor
from pathfinder.process.cleaners import SurveyCleaner
from pathfinder.process.trackers import SurveyTracker
from pathfinder.process.data import SurveyData
#######################
# Results base object #
#######################


class Result:

    def __init__(
            self,
            path: str or Path = None,
            file: str or Path = None
    ):

        self.path = path
        self.file = file

        self.data = None

    def __add__(self, other):

        # Enable sum() with __radd__
        # by skipping starting condition:
        if other == 0:
            return self
        else:
            if not isinstance(other, Result):
                raise TypeError('Only other Result objects can be added.')
            self.path = None
            self.file = None
            self.data += other.data

            return self

    def __radd__(self, other):

        return self.__add__(other)

    def complete(self):
        """ Overwritten in subclasses """
        pass

################################
# Result objects for pipelines #
################################


class SurveyResult(Result, SurveyProcessor, SurveyCleaner):

    """ Class for cleaning and storing results from pf-core/pf-survey """
    def __init__(
            self,
            path: str or Path = None,
            file: str or Path = None
    ):
        """ Public access to cleaning parameters on initialisation """

        Result.__init__(self, path=path, file=file)

        self.data = SurveyData()
        self.tracker = SurveyTracker()

        if file:
            self.tracker.read(
                file=Path(file) if isinstance(file, str) else file
            )

        SurveyProcessor.__init__(
            self, path=self.path, data=self.data, tracker=self.tracker
        )   # Exposes parsed and processed self.data

        SurveyCleaner.__init__(
            self, data=self.data, tracker=self.tracker
        )   # Uses self.data to subset

    def complete(self):

        """ Complete simply returns a subset of the data according
        to a completion criterion, for instance here this returns
        the parsed data subset by the IDs that passed MLST and Kraken
        """

        df = self.tracker.as_dataframe()

        # Complete = intersection of MLST and Kraken

        iids_mlst = set(
            df.loc[df['mlst'], 'mlst'].index
        )
        iids_kraken = set(
            df.loc[df['kraken'], 'kraken'].index
        )
        complete = list(
            iids_mlst.intersection(iids_kraken)
        )

        self.data = self.data.subset(iids=complete)

        for attr, df in self.data:
            self.tracker.complete += [df.index.unique().tolist()]

        print(f'Total number of SurveyResult-complete genomes: {len(complete)}')

    def print_mlst_summary(self, top=100, min_count=100, plot: str = None):

        pandas.set_option('display.max_rows', top)

        counts = self.data.mlst['sequence_type'].value_counts()

        if min_count:
            counts = counts[counts >= min_count]

        print(
            counts.nlargest(top)
        )

    def plot_mlst_summary(self, top=10):

        self.data.mlst['sequence_type'].value_counts()

