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

    def update_iid(self):

        """ Update list of unique IIDs in current data container """

        iids = pandas.Series(
            [iid for attr, df in self for iid in df.index.unique().tolist()]
        ).unique()

        self.iid = pandas.DataFrame(data={'iid': iids}, index=iids)

    def remove(self, process_remove_dict):

        for process, df in self:

            # TODO: Patched, needs proper thought.
            try:
                process_remove_dict[process]
            except KeyError:
                # Skip self.iid, self.fasta, self.fastq
                continue

            removed_df = df.drop(
                process_remove_dict[process], errors='ignore'
            )
            setattr(self, process, removed_df)

        self.update_iid()

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
        for this ID is introduced in to the data frame.
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

            print(largest)
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

        for fpath in fasta_dir.glob(f'*{extension}'):
            iid = fpath.name.replace(extension, '')
            if iid in self.fasta.index:
                self.fasta.at[iid, 'fasta'] = fpath.resolve()
            else:
                pass
                # Log failure if IID not in DF index:
                # print(f'{iid} not in fasta index.')

    def _add_fastq(self, fastq_dir: Path, forward_tail: str = '_1',
                   reverse_tail: str = '_2', extension: str = '.fq.gz'):

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

        Path(result_path) if isinstance(result_path, str) else result_path

        if fasta_dir:
            self._add_fasta(result_path / fasta_dir)
        if fastq_dir:
            self._add_fastq(result_path / fastq_dir)


@dataclass
class SurveyData(ResultData):

    mlst: pandas.DataFrame = pandas.DataFrame()
    kraken: pandas.DataFrame = pandas.DataFrame()
    mash: pandas.DataFrame = pandas.DataFrame()

    abricate_resistance: pandas.DataFrame = pandas.DataFrame()
    abricate_virulence: pandas.DataFrame = pandas.DataFrame()
    abricate_plasmid: pandas.DataFrame = pandas.DataFrame()

    mykrobe_phenotype: pandas.DataFrame = pandas.DataFrame()
    mykrobe_genotype: pandas.DataFrame = pandas.DataFrame()

    def __copy__(self):

        sd = SurveyData()
        for attr, data in self:
            setattr(sd, attr, data)
        return sd

    def get_file_paths(
            self,
            result_path: str or Path,
            fasta_dir: str = 'skesa',
            fastq_dir: str = 'trimmomatic'
    ):
        super(SurveyData, self).get_file_paths(
            result_path, fasta_dir, fastq_dir
        )

