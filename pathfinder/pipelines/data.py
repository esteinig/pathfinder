import pandas
from pandas.core.groupby import DataFrameGroupBy
from numpy import nan
from pathlib import Path
from dataclasses import dataclass
import copy
import shutil
from tqdm import tqdm


#################################
# Data containers for Nextflows #
#################################


@dataclass
class ResultData:
    """ Methods for data containers. """

    iid: pandas.DataFrame = pandas.DataFrame(
        columns=['iid']
    )
    fasta: pandas.DataFrame = pandas.DataFrame(
        columns=['fasta']
    )
    fastq: pandas.DataFrame = pandas.DataFrame(
        columns=['forward', 'reverse']
    )

    def __copy__(self):
        """ Overwritten in subclasses """

    def __iter__(self):

        for attr, df in self.__dict__.items():
            yield attr, df

    def __add__(self, other):

        for attr, data in self:
            data_other = getattr(other, attr)

            merged = pandas.concat(
                [data, data_other], sort=False
            )

            setattr(self, attr, merged)

        self.update_iid()

        return self

    def __getitem__(self, iid):

        return self.iid[iid]

    def update_iid(self, iids=None):

        """ Update list of unique IIDs in current data container """

        if not iids:
            iids = pandas.Series(
                [iid for attr, df in self for iid in df.index.unique().tolist()]
            ).unique()

        self.iid = pandas.DataFrame(data={'iid': iids}, index=iids)

    def remove(self, remove_dict: dict, retain: bool = False):

        for process, df in self:

            # TODO: Patched, needs proper thought.
            try:
                remove_dict[process]
            except KeyError:
                # Skip self.iid, self.fasta, self.fastq
                continue

            if retain:
                removed_df = df[
                    df.index.isin(
                        remove_dict[process]
                    )
                ]
            else:
                removed_df = df.drop(
                    remove_dict[process], errors='ignore'
                )

            setattr(self, process, removed_df)

        self.update_iid()

    def complete(
            self,
            at: list = ('kraken', 'mlst')
    ):

        """ Return data object with the set intersection between all IIDs
        in the analyses results; usually after filtering, can be narrowed
        to only some results by default to Kraken and MLST """

        if at:
            allowed = at
        else:
            allowed = [process for process, _ in self]

        iids = [
            set(df.index.tolist()) for process, df in self
            if len(df) > 0 and process in allowed
        ]

        for li in iids:
            print('Len IID data:', len(li))

        intersect = list(set.intersection(*iids))

        print('Intersect, keep this many:', len(intersect))

        self.remove(
            {process: intersect for process, _ in self},
            retain=True
        )

    def write(self, outdir: str or Path) -> None:

        outdir = Path(outdir).resolve() \
            if isinstance(outdir, str) else outdir.resolve()

        outdir.mkdir(parents=True, exist_ok=True)

        for process, data in self:
            data.to_csv(
                (outdir / f'{process}.csv'),
                index=True, header=True, sep=','
            )

    def read(self, results: str or Path):

        results = Path(results).resolve() \
            if isinstance(results, str) else results.resolve()

        for file in results.glob('*.csv'):
            setattr(self, file.stem, pandas.read_csv(
                file, index_col=0, header=0, sep=','
            ))

        self.update_iid()

    def subset(self, iids, inplace=False):
        """ Subset the Data by a list of isolate IDs.

        If IDs are not in the data frame index, a row of null values
        for this ID is introduced into the data frame.
        """

        if inplace:
            sub = self
        else:
            sub = copy.copy(self)

        for attr, df in self:
            subset = df[df.index.isin(iids)]
            missing = list(set(iids).difference(set(subset.index)))
            if missing:
                subset = pandas.concat(
                    (
                        subset,
                        pandas.DataFrame(index=missing, columns=df.columns)
                     )
                )

            setattr(sub, attr, subset)

        sub.update_iid()

        return sub

    def groupby(self, attr, column, set_index=False, **kwargs) -> [DataFrameGroupBy]:

        """ Map the column values to each corresponding index
        in all other data attributes and group the data frames
        by the column values. This allows e.g. to group all data
        sequence type etc.
        """

        if not hasattr(self, attr):
            raise AttributeError(
                f'ResultData has no attribute: {attr}'
            )

        # Map sequence type to index in each DataFrame as new column
        # and then group every attribute data frame by this column:

        data = getattr(self, attr)
        data_series = data[column]

        for attr, df in self:
            df[column] = df.index.map(data_series)

            if set_index:
                df = df.set_index([df.index, column])
                setattr(self, attr, df)

            groups = df.groupby(by=column, **kwargs)
            yield groups, attr

    def select(self, attr, column, min_count=None,
               sample=None, values=None):

        data = getattr(self, attr)

        # TODO: Both

        if values:
            iids = data[
                data[column].isin(values)
            ].index.unique().tolist()

        elif min_count:
            counts = data[column].value_counts()
            largest = counts[counts > min_count]

            subset = largest.index.to_series()  # Index: unique values counted

            sub = data[
                data[column].isin(subset)
            ]

            iids = sub.index.unique().tolist()

            if sample:

                iids = pandas.concat([
                    group.sample(n=sample) for i, group in sub.groupby(column)
                ]).index

        else:
            iids = data[column].index.unique().tolist()

        return self.subset(iids)

    def link_fasta(self, fdir: str or Path, symlink: bool = True,
                   progbar: bool = True, index: dict = None):

        fdir = self._check_dir(fdir)

        # Can we predict read length distribution from a smear
        # of gel run of HMW DNA

        for fasta in tqdm(
                self.fasta.fasta,
                disable=not progbar,
                total=len(self.fasta)
        ):
            if index:
                try:
                    name = index[Path(fasta).stem] + '.fasta'
                except TypeError:
                    print(index[Path(fasta).stem])
                    continue
                except KeyError:
                    print('Could not find entry in FASTA Index.')
                    continue
            else:
                name = Path(fasta).name

            if symlink:
                (fdir / name).symlink_to(fasta)
            else:
                shutil.copy(
                    fasta, str(fdir / name)
                )

    def link_fastq(self, fdir: str or Path, symlink: bool = True,
                   progbar: bool = True):

        fdir = self._check_dir(fdir)

        for fwd, rv in tqdm(
            self.fastq.itertuples(index=False),
            disable=not progbar,
            total=len(self.fastq)
        ):
            for fastq in (fwd, rv):
                if symlink:
                    (fdir / Path(fastq).name).symlink_to(fastq)
                else:
                    shutil.copy(
                        fastq, str(fdir)
                    )

    @staticmethod
    def _check_dir(fdir: str or Path) -> Path:

        fdir = Path(fdir) if isinstance(fdir, str) else fdir

        if not fdir.exists():
            fdir.mkdir(parents=True)

        return fdir

    # File operations

    def _add_fasta(self, fasta_dir: Path, extension: str = '.fasta'):

        """ Add file paths to data frame in attribute: `fasta` """

        # Clear FASTA, weirdly in loops files accumulate
        # even if new instance of this is created

        self.fasta = pandas.DataFrame()

        for fpath in fasta_dir.glob(f'*{extension}'):
            iid = fpath.name.replace(extension, '')
            self.fasta.at[iid, 'fasta'] = fpath.resolve()

    def _add_fastq(self, fastq_dir: Path, forward_tail: str = '_1',
                   reverse_tail: str = '_2', extension: str = '.fq.gz'):
        

        self.fastq = pandas.DataFrame()

        for fpath in fastq_dir.glob(f'*{forward_tail}{extension}'):
            iid = fpath.name.replace(f'{forward_tail}{extension}', '')
            if iid in self.fastq.index:
                self.fastq.at[iid, 'forward'] = fpath.resolve()
            else:
                pass

        for fpath in fastq_dir.glob(f'*{reverse_tail}{extension}'):
            iid = fpath.name.replace(f'{reverse_tail}{extension}', '')
            if iid in self.fastq.index:
                self.fastq.at[iid, 'reverse'] = fpath.resolve()
            else:
                pass

    def get_file_paths(
        self,
        result_path: str or Path,
        fasta_dir: str = None,
        fastq_dir: str = None
    ) -> None:

        if fasta_dir is None and fastq_dir is None:
            raise ValueError(
                'Please specify a directory name for FASTA or FASTQ.'
            )

        print('Getting file paths:', result_path)

        Path(result_path) if isinstance(result_path, str) else result_path

        if fasta_dir:
            self._add_fasta(result_path / fasta_dir)
        if fastq_dir:
            self._add_fastq(result_path / fastq_dir)

        self.update_iid()

    def by_iid(self) -> pandas.DataFrame:

        """ Returns a bool summary of IIDs in each process """

        exclude_attrs = ('fasta', 'fastq', 'iid')

        df = pandas.DataFrame(
            index=self.iid.iid,
            columns=[attr for attr, _ in self if attr not in exclude_attrs]
        )

        for attr, data in self:
            if attr in exclude_attrs:
                continue
            for iid in self.iid.iid:
                if iid in data.index:
                    df.at[iid, attr] = True
                else:
                    df.at[iid, attr] = False

        return df


@dataclass
class SurveyData(ResultData):

    mlst: pandas.DataFrame = pandas.DataFrame()
    kraken: pandas.DataFrame = pandas.DataFrame()
    mash: pandas.DataFrame = pandas.DataFrame()
    kleborate: pandas.DataFrame = pandas.DataFrame()

    abricate_resistance: pandas.DataFrame = pandas.DataFrame()
    abricate_virulence: pandas.DataFrame = pandas.DataFrame()
    abricate_plasmid: pandas.DataFrame = pandas.DataFrame()

    mykrobe_phenotype: pandas.DataFrame = pandas.DataFrame()
    mykrobe_genotype: pandas.DataFrame = pandas.DataFrame()
    mykrobe_lineage: pandas.DataFrame = pandas.DataFrame()

    def __copy__(self):

        sd = SurveyData()
        for attr, data in self:
            setattr(sd, attr, data)
        return sd

    @property
    def empty(self):

        check = [data.empty for attr, data in self]
        if all(check):
            return True
        else:
            return False

    def sketchy(
        self,
        config: dict = None,
        lineage: dict = None,
        genotype: dict = None,
        susceptibility: dict = None,
        sep1: str = '-',
        sep2: str = '',
    ) -> pandas.DataFrame:

        """ Create a data frame for MinHash sketching in Sketchy """

        if config:
            try:
                lineage = config['lineage']
            except KeyError:
                lineage = None

            try:
                genotype = config['genotype']
            except KeyError:
                genotype = None

            try:
                susceptibility = config['susceptibility']
            except KeyError:
                susceptibility = None

        data = dict(
            lineage=self._get_sketchy_data(lineage, sep=sep1),
            genotype=self._get_sketchy_data(genotype, sep=sep1),
            susceptibility=self._get_sketchy_data(susceptibility, sep=sep2),
            fasta=self.fasta.fasta
        )

        # Cleaning step if using lineage from Mykrobe, arcane #TODO
        lineage = data['lineage']
        data['lineage'] = lineage.loc[~lineage.index.duplicated()]

        df = pandas.DataFrame().from_dict(data)

        return df

    def _get_sketchy_data(self, data: dict or None, sep: str):

        """
        Extract and concatenate data across result data frames for Sketchy
        """

        if data is None:
            return pandas.Series(index=self.iid.iid)

        dfs = []
        for attr, columns in data.items():
            df = getattr(self, attr)

            if isinstance(columns, str):
                columns = [columns]

            if isinstance(columns, pandas.Index or pandas.Series):
                columns = df.columns.tolist()

            # Get all:
            if len(columns) == 0:
                columns = df.columns.tolist()

            test = df[columns].apply(
                lambda row: f'{sep}'.join(
                    row.dropna().values.astype(str)
                ), axis=1
            )
            d = pandas.DataFrame({attr: test})

            # print(d.index[:5], df.index[:5])
            # print(len(d), len(df))

            dfs += [d]

        complete = pandas.concat(dfs, axis=1)

        return complete.apply(lambda row: f'{sep}'.join(
            row.dropna().values.astype(str)
        ), axis=1)




