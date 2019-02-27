import pandas
from colorama import Fore
from textwrap import dedent
from pathfinder.process.settings import SurveyCleanerSetting
from pathfinder.process.data import SurveyData
from pathfinder.process.trackers import SurveyTracker
from numpy import unique
from pathlib import Path


class Cleaner:
    """ Functions to clean Parser + Processor results """

    def __init__(self, data, settings=None):

        self.data = data
        self.settings = settings  # Settings Dataclass types overwritten

        self.remove = dict()
        self.remove_flattened = list()

        self.tracker = None

    def clean(self):
        """ Overwritten in subclasses """

        pass

    # Main cleaning functions for mix-and-matching subclasses

    def clean_kraken(
            self, data: pandas.DataFrame
    ) -> list:

        """ Clean and condense the taxonomic read assignments by Kraken2

        Inspects each isolate's results from the agglomerated data frame
        and identify the following:

            .. :py:attr:`pathfinder.results.SurveySettings.top`

                Number of species to summarize by largest percentage of
                assigned reads.

            .. py:attr:`pathfinder.results.SurveySettings.contamination`

                Threshold for removing contamination with this percentage of
                reads classified in top species besides the most represented
                species.

            .. py:attr:`pathfinder.results.SurveySettings.purity`

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

                if g['percent'].iloc[0] < self.settings.purity*100:
                    to_remove.append(iid)
                    uncertain += 1

                if (g['percent'][1:] > self.settings.contamination*100).any():
                    contaminated += 1
                    to_remove.append(iid)

                if self.settings.species:
                    if isinstance(self.settings.species, str):
                        column = 'taxonomy'
                    else:
                        column = 'taxid'

                    if g[column].iloc[0] != self.settings.species:
                        to_remove.append(iid)
                        misidentified += 1

        return list(
            pandas.Series(to_remove).unique()
        )

    def clean_mykrobe(self, data: pandas.DataFrame) -> list:

        return list()

    def clean_mash(self, data: pandas.DataFrame) -> list:

        return list()

    def clean_abricate(self, data: pandas.DataFrame) -> list:

        """ Subset and reduce Abricate data frames """

        return list()

    def clean_mlst(self, data: pandas.DataFrame) -> list:

        if self.settings.complete_lineage:
            return data[data['sequence_type'] == '-'].index.tolist()
        else:
            return list()


class SurveyCleaner(Cleaner):

    def __init__(self, data: SurveyData, tracker: SurveyTracker = None):

        Cleaner.__init__(self, data=data)

        self.tracker = tracker

    def clean(
            self,
            settings: SurveyCleanerSetting = None, **kwargs
    ):
        """ Explicit cleaning subclass for pf-core/pf-survey Nextflow """

        if settings is None:
            self.settings = SurveyCleanerSetting(**kwargs)
        else:
            self.settings = settings

        self.remove = {
            'kraken': self.clean_kraken(
                self.data.kraken
            ),
            'mlst': self.clean_mlst(
                self.data.mlst
            ),
            'mash': self.clean_mash(
                self.data.mash
            ),
            'abricate_resistance': self.clean_abricate(
                self.data.abricate_resistance
            ),
            'abricate_virulence': self.clean_abricate(
                self.data.abricate_virulence
            ),
            'abricate_plasmid': self.clean_abricate(
                self.data.abricate_plasmid
            ),
            'mykrobe_genotype': self.clean_mykrobe(
                self.data.mykrobe_genotype
            ),
            'mykrobe_phenotype': self.clean_mykrobe(
                self.data.mykrobe_phenotype
            ),
        }

        self.remove_flattened = pandas.Series(
            [iid for _, cleaned_iids in self.remove.items()
             for iid in cleaned_iids]
        ).unique()

        print(
            f'Cleaned {Fore.YELLOW}{len(self.remove_flattened)}'
            f'{Fore.RESET} genomes.'
        )

        self.data.remove(self.remove)

        self._track_cleaned()

    def _track_cleaned(self):

        for process in self.tracker.process:
            data = getattr(self.data, process)
            self.tracker.cleaned.append(
                data.index.unique().tolist()
            )

