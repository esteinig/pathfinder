"""

Parsing and process of data from Nextflow follows:

    parse - aggregate - process

That is, output files from the pipelines are first parsed into
:py:class:`pandas.DataFrame`, then aggregated into a larger
:py:class:`pandas.DataFrame` of all the samples processed in
the workflow, and then data may be subset or processed to make
sense of it.

This structure is handled by the following superclasses:

    .. py:class:: pathfinder.results.Processor:

        Super class contains anonymous `parse - process - agglomerate`
        methods for running subclass process of pipeline data.

        Inherited by :py:class:`nanopath.results.Result`.

    .. py:class:: pathfinder.results.Result:

        Super class for generic results process functions. I

Superclasses are used nd the following pipeline-specific subclasses:

    .. py:class:: pathfinder.results.SurveyProcessor:

        Handles the

        Inherited by :py:class:`nanopath.results.SurveyProcessor`.

    .. py:class:: pathfinder.results.SurveyResult:

        Inherits from :py:class:`nanopath.results.SurveyProcessor`.

"""

import os
import dataclasses
import tqdm
import pandas
from colorama import Fore

from pathlib import Path
from pathfinder.process.parsers import Parser
from pathfinder.process.data import SurveyData
from pathfinder.process.trackers import SurveyTracker
from pathfinder.process.settings import SurveyProcessSetting

Y = Fore.YELLOW
R = Fore.RESET

###################################
# Processor classes for Nextflows #
###################################


class Processor:

    """ Super class for process any result from Nextflow"""

    def __init__(
            self,
            path: str or Path = None,
            file: str or Path = None,
    ):
        """ Init with list of IDs from Nextflow output

        This collects all process functions, similar to Parser and Cleaner
        classes - the reason being that users might want to access any of
        the parse - process - clean methods for each piece of software and build
        their own class that parses - settings - cleans pipeline output. By
        bundling all base methods into these objects, pipeline specific result
        classes only need to inherit from Parser, Processors and Cleaners.
        """

        self.path = Path(path) if isinstance(path, str) else path
        self.file = Path(file) if isinstance(file, str) else file

        self.data = None  # Assigned in subclasses, subclass specific Dataclass
        self.tracker = None  # Assigned in subclasses, subclass specific Tracker

    def get_result(self) -> []:
        """ Overwritten with process methods in subclasses """
        return []

    def process(
            self,
            progbar: tqdm.tqdm = None,
    ) -> None:
        """ Explicit parsing and process results from survey Nextflow """

        for process, df in self.get_result():
            setattr(self.data, process, df)
            self.tracker.process += [process]

            if progbar:
                progbar.update()

        self.data.update_iid()

    def _parse_and_process(
            self,
            files: [],
            parse_func=None,
            process_func=None,
            aggregate: bool = True,
            remove: str or list = None,
            subset: str = None, *args, **kwargs
    ) -> dict or pandas.DataFrame:

        """Anonymous function to parse results based on a list of files,
        the name components to be stripped for extracting IDs, and the
        parse method of this class to return a dictionary of:

            .. key

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
        :param aggregate:
            aggregate all results into a single dataframe

        :returns parsed results (dict) or aggregated results (dataframe)
        """

        files = list(files)

        if subset:
            files = retain_files(files, retain=subset)

        self.tracker.files.append(files)

        # Parse results from many files by ID:
        data = {
            get_id_from_fname(file, remove=remove): parse_func(file, *args)
            for file in files
        }

        self.tracker.parsed.append(
            list(data.keys())
        )

        # TODO: always aggregate!
        # Aggregate the data into DataFrames:
        data = self._aggregate(data)

        self.tracker.aggregated.append(
            data.index.unique().tolist()
        )

        if process_func:
            # Process the DataFrames into reduced formats
            data = process_func(data, **kwargs)
            self.tracker.processed.append(
                data.index.unique().tolist()
            )
        else:
            self.tracker.processed.append(None)

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

    #######################################
    # Process methods for SurveyProcessor #
    #######################################

    @staticmethod
    def process_fastqc(df):
        """Process and clean FastQC dataframes"""

        df = df[df.category != "Per tile sequence quality"]
        id_grouped = df.groupby(df.index)
        id_groups = [
            id_grouped.get_group(id_group) for id_group in id_grouped.groups
        ]

        dataframes = []
        for df in id_groups:
            read_grouped = df.groupby('file')

            forward, reverse = [
                read_grouped.get_group(read) for read in read_grouped.groups
            ]

            frc = forward.reset_index() \
                .merge(reverse, on="category", how='outer') \
                .set_index('ID') \
                .drop(["file_x", "file_y"], axis=1) \
                .rename(
                    {"alert_x": "forward", "alert_y": "reverse"}, axis=1
                )

            # Pivoting to get the data row-wise:
            forward = frc.pivot(
                index=None, columns='category', values=["forward"]
            )
            forward.columns = forward.columns.droplevel(0)
            reverse = frc.pivot(
                index=None, columns='category', values=["reverse"]
            )
            reverse.columns = reverse.columns.droplevel(0)

            frc = pandas.concat([forward, reverse]).assign(
                reads=["forward", "reverse"]
            )

            dataframes.append(frc)

        return pandas.concat(dataframes)

    @staticmethod
    def process_abricate(df):
        # Remove alleles from gene names:
        df.gene = df.gene.str.rsplit('_').str[0]
        return df

    @staticmethod
    def process_kraken(df):
        df.taxonomy = df.taxonomy.str.strip()
        return df


class SurveyProcessor(Parser, Processor):

    def __init__(
            self,
            path: Path,
            data: SurveyData = None,
            tracker: SurveyTracker = None
    ):

        Parser.__init__(self)
        Processor.__init__(
            self, path=path
        )

        ##########################
        # Initialisation Process #
        ##########################

        self.data = data
        self.tracker = tracker

        self.processes = {
            'kraken': SurveyProcessSetting(
                files=list(
                    (self.path / 'kraken').glob('*.report')
                ),
                parse_func=self.parse_kraken,
                process_func=self.process_kraken,
                remove='.report',
                aggregate=True
            ),
            'mlst': SurveyProcessSetting(
                files=list(
                    (self.path / 'mlst').glob('*.tab')
                ),
                parse_func=self.parse_mlst,
                process_func=None,
                remove='.tab',
                aggregate=True
            ),
            'mash': SurveyProcessSetting(
                files=list(
                    (self.path / 'mash').glob('*.mash.tab')
                ),
                parse_func=self.parse_mash,
                process_func=None,
                remove='.mash.tab',
                aggregate=True
            ),
            'mykrobe_phenotype': SurveyProcessSetting(
                files=list(
                    (self.path / 'mykrobe').glob('*.json')
                ),
                parse_func=self.parse_mykrobe_susceptibility,
                process_func=None,
                remove='.json',
                aggregate=True
            ),
            'mykrobe_genotype': SurveyProcessSetting(
                files=list(
                    (self.path / 'mykrobe').glob('*.json')
                ),
                parse_func=self.parse_mykrobe_genotype,
                process_func=None,
                remove='.json',
                aggregate=True
            ),
            'abricate_resistance': SurveyProcessSetting(
                files=list(
                    (self.path / 'abricate' / 'resfinder').glob('*.tab')
                ),
                parse_func=self.parse_abricate,
                process_func=self.process_abricate,
                remove='.tab',
                aggregate=True
            ),
            'abricate_virulence': SurveyProcessSetting(
                files=list(
                    (self.path / 'abricate' / 'vfdb').glob('*.tab')
                ),
                parse_func=self.parse_abricate,
                process_func=self.process_abricate,
                remove='.tab',
                aggregate=True
            ),
            'abricate_plasmid': SurveyProcessSetting(
                files=list(
                    (self.path / 'abricate' / 'plasmidfinder').glob('*.tab')
                ),
                parse_func=self.parse_abricate,
                process_func=self.process_abricate,
                remove='.tab',
                aggregate=True
            )
        }

    def get_result(self) -> (str, pandas.DataFrame):

        for process, process_setting in self.processes.items():
            print(f'Processing: {Y}{process}{R}')
            yield process, self._parse_and_process(
                **dataclasses.asdict(process_setting)
            )


###########################
# Module helper functions #
###########################


def get_content(path):
    """Get content by directory / files if it contains files """
    content = {}
    for directory, subdirs, files in os.walk(path):
        if files:
            content[Path(directory)] = [
                Path(directory) / Path(file) for file in files
            ]

    return content


def get_id_from_fname(fname: str or Path, remove: str or list = None):
    """Helper function to deconstruct a filename
    since there is no scheme to IDs.
    :param fname: file basename.
    :param remove: substrings to be removed."""

    if not remove:
        return fname

    if isinstance(fname, Path):
        fname = str(fname.name)
    if isinstance(remove, list):
        for r in remove:
            fname = fname.replace(r, '')
        return fname
    elif isinstance(remove, str):
        return fname.replace(remove, '')
    else:
        raise ValueError


def retain_files(files: list, retain: str) -> list:
    """Helper function to retain files if sub str
    is contained in file name
    :param files: list of file paths
    :param retain: str in file name, to retain file
    """
    retained = []
    for file in files:
        if isinstance(file, Path):
            fname = str(file.name)
        else:
            fname = os.path.basename(file)
        if retain in fname:
            retained.append(file)

    return retained

