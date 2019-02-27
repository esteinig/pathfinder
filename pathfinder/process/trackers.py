import pandas

from pathlib import Path
from dataclasses import dataclass
from colorama import Fore

Y = Fore.YELLOW
R = Fore.RESET


@dataclass
class TrackerData:
    """ Dataclass for stages of tracking isolate IDs """
    def __init__(self):

        self.iid: list = list()
        self.process: list = list()
        self.files: [list] = list()
        self.parsed: [list] = list()
        self.aggregated: [list] = list()
        self.processed: [list] = list()
        self.cleaned: [list] = list()
        self.complete: [list] = list()

    def as_dataframe(self):

        df = pandas.DataFrame(
            index=self.iid,
            columns=self.process
        )

        for i, col in enumerate(self.process):
            for iid in self.iid:
                if iid in getattr(self, 'cleaned')[i]:
                    df.at[iid, col] = True
                else:
                    df.at[iid, col] = False

        return df

    def read(
            self,
            file: Path,
            sep: str = ',',
            iid: str = 'index',
            header: int = 0
    ) -> pandas.DataFrame:

        if iid == 'index':
            df = pandas.read_csv(file, sep=sep, header=header, index_col=1)
        else:
            df = pandas.read_csv(file, sep=sep, index_col=iid, header=header)

        self.iid = df.index.tolist()

        return df


@dataclass
class SurveyTracker(TrackerData):

    def __init__(self):

        TrackerData.__init__(self)

    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value

    @staticmethod
    def _print(rows, fmt):

        print('\n', TablePrinter(fmt, ul='=')(rows), '\n')

    def print_steps(self):

        fmt = [
            ('Process', 'process', 20),
            ('Files', 'files', 15),
            ('Parsed IDs', 'parsed', 15),
            ('Aggregated IDs', 'aggregated', 15),
            ('Processed IDs', 'processed', 15),
            ('Cleaned IDs', 'cleaned', 15),
            ('Complete IDs', 'complete', 15)
        ]

        rows = []
        for i, process in enumerate(self.process):
            data = {'process': process}
            for col in fmt[1:]:
                attr_name = col[1]
                value = getattr(self, attr_name)
                try:
                    if value[i] is not None:
                        data[attr_name] = len(value[i])
                    else:
                        data[attr_name] = 'null'

                except IndexError:
                    data[attr_name] = 'null'

            rows.append(data)

        self._print(rows=rows, fmt=fmt)


class TablePrinter(object):
    """ Print a list of dicts as a table """

    def __init__(self, fmt, sep=' ', ul=None):
        """
        @param fmt: list of tuple(heading, key, width)
                    heading: str, column label
                    key: dictionary key to value to print
                    width: int, column width in chars
        @param sep: string, separation between columns
        @param ul: string, character to underline column label,
         or None for no underlining
        """

        super(TablePrinter, self).__init__()

        self.fmt = str(sep).join(
            '{lb}{0}:^{1}{rb}'.format(key, width, lb='{', rb='}')
            for heading, key, width in fmt
        )

        self.head = {key: heading for heading, key, width in fmt}

        self.ul = {
            key: str(ul)*width for heading, key, width in fmt
        } if ul else None

        self.width = {key: width for heading, key, width in fmt}

    def row(self, data):
        return self.fmt.format(
            **{k: str(data.get(k, ''))[:w] for k, w in self.width.items()})

    def __call__(self, data_list):
        _r = self.row
        res = [_r(data) for data in data_list]
        res.insert(0, _r(self.head))
        if self.ul:
            res.insert(1, _r(self.ul))
        return '\n'.join(res)
