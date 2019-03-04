import time
import uuid
import pandas
import urllib.request

from pathfinder.utils import stamp, get_genome_sizes
from pathfinder.exceptions import SurveyError

import shlex
import subprocess

from io import StringIO

from tqdm import tqdm
from pathlib import Path
from pandas.errors import EmptyDataError


class MiniAspera:

    def __init__(self, force=False):

        self.port = 33001
        self.limit = 1024

        self.force = force

        self.ascp = Path.home() / ".aspera/connect/bin/ascp"
        self.key = Path.home() / ".aspera/connect/etc/asperaweb_id_dsa.openssh"

        self.fasp = "era-fasp@fasp.sra.ebi.ac.uk:"

    def download_batch(
            self,
            file,
            outdir: str = ".",
            limit_download: int = None
    ):

        batch = self.read_batch(file)

        outpath = Path(outdir)
        if not outpath.is_dir():
            outpath.mkdir(parents=True, exist_ok=True)

        if limit_download:
            batch = batch.iloc[0:limit_download]

        with tqdm(total=len(batch)) as pbar:
            pbar.set_description("Downloading batch")

            for i, fastq in batch.iterrows():
                fq1_address = fastq["ftp_1"].replace(
                    "ftp.sra.ebi.ac.uk", self.fasp
                )
                fq1_path = Path(outdir) / Path(fq1_address).name
                self.download(
                    address=fq1_address, outfile=fq1_path, force=self.force
                )

                if fastq["ftp_2"]:
                    fq2_address = fastq["ftp_2"].replace(
                        "ftp.sra.ebi.ac.uk", self.fasp
                    )
                    fq2_path = Path(outdir) / Path(fq2_address).name
                    self.download(
                        address=fq2_address, outfile=fq2_path, force=self.force
                    )

                pbar.update(1)

    def download(self, address, outfile, force=False):

        # Skip existing files
        if not force and outfile.exists():
            print(f"File exists: {outfile}")
            return

        cmd = shlex.split(
            f"{self.ascp} -QT -l {str(self.limit)}m"
            f" -P{str(self.port)} -i {self.key} -q "
            f"{address} {outfile}"
        )

        try:
            subprocess.call(cmd)
        except subprocess.CalledProcessError:
            print("Error in subprocess.")
            raise  # handle errors in the called executable
        except OSError:
            print("Executable not found.")
            raise  # executable not found

    @staticmethod
    def read_batch(file):

        return pandas.read_csv(file)


class Survey:

    """
    Simple accessor class wrapping wget to pull short-read data from the ENA.
    """

    def __init__(self, outdir=None):

        self.outdir = outdir

        self.url_result = "read_run"
        self.url_display = "report"
        self.url_query = "https://www.ebi.ac.uk/ena/data/" \
                         "warehouse/search?query="
        self.url_fields = "run_accession,tax_id,fastq_ftp,fastq_bytes," \
                          "read_count,base_count," \
                          "instrument_platform,instrument_model," \
                          "library_layout,library_source," \
                          "library_strategy,sample_accession,study_accession"

        # TODO add study accession

        self.results = dict()
        self.query = pandas.DataFrame()

    def query_to_csv(self, csv_file="query.csv", query_results=None):

        if not query_results:
            self.query.to_csv(csv_file)
        else:
            query_results.to_csv(csv_file)

    def query_from_csv(self, file):

        self.query = pandas.read_csv(file)

        return self.query

    def parse_biosample(self):

        pass

    def display(self):

        pass

    def filter_query(self, ops: str):
        """ Filter the query dataframe by
        applying pandas query (boolean
        filter operations) string """
        self.query = self.query.query(ops)

    @staticmethod
    def batch_output(batches, outdir="batches", exist_ok=True):

        outdir = Path.cwd() / outdir
        for i, batch in enumerate(batches):
            batch_dir = outdir / f"batch_{i}"
            batch_csv = batch_dir / f"batch_{i}.csv"

            batch_dir.mkdir(parents=True, exist_ok=exist_ok)
            batch.to_csv(batch_csv)

            yield batch_dir, batch_csv

    def batch(self, query=None, batch_size=None, max_gb=None):

        if query is None:
            query = self.query

        if max_gb:
            running_size = []
            batches = []
            start = 0
            for i, entry in query.iterrows():
                running_size.append(float(entry["size"])/1000)
                if sum(running_size) > max_gb:
                    batches += [query[start:i]]
                    start = i
                    running_size = []
            return batches
        elif batch_size:
            return [
                query[i:i + batch_size] for i in range(
                    0, query.shape[0], batch_size
                )
            ]
        else:
            raise ValueError(
                "Either maximum gigabytes or batch size must be set."
            )

    def query_ena(self, species="Staphylococcus aureus", scheme="illumina",
                  study: str = None, term: str = None) -> (dict, str):

        """
        Search the ENA warehouse for raw sequence reads with the
        following default parameters to conform to current implementations
        of WGS analysis pipelines.
        """

        # Format queries correctly:
        scheme = scheme.lower()

        # At the moment, restrict queries to schemes:
        if scheme == "illumina":

            platform = "ILLUMINA"
            source = "GENOMIC"
            layout = "PAIRED"
            strategy = "WGS"

        elif scheme == "nanopore":

            platform = "OXFORD_NANOPORE"
            source = "GENOMIC"
            layout = "SINGLE"
            strategy = "WGS"

        else:
            raise SurveyError(
                "Require pre-defined WGS scheme: 'Illumina' or 'Nanopore'"
            )

        if species:
            term = self._construct_species_query(
                species, platform, source, layout, strategy
            )
        elif study:
            term = self._construct_study_query(study)
        elif term:
            pass
        else:
            raise ValueError(
                "Need to specify either species, study accession, "
                "or custom search term for Survey.")

        url = f"{self.url_query}{term}&result={self.url_result}" \
              f"&fields={self.url_fields}&" \
            f"display={self.url_display}".replace(" ", "%20")

        query_results = StringIO(urllib.request.urlopen(url)
                                 .read().decode('utf-8'))

        try:
            df = pandas.read_csv(query_results, sep="\t")
        except EmptyDataError:
            raise SurveyError(f"No results were returned for query: {url}")

        query_results = self._sanitize_ena_query(df, url)

        self.results[time.time()] = query_results
        self.query = query_results

        return query_results, term

    @staticmethod
    def _sanitize_ena_query(df, url) -> pandas.DataFrame:

        # Drop rows with missing FTP links:
        df = df.dropna(subset=["fastq_ftp"])

        if df.empty:
            raise SurveyError(f"Query results are empty, check your "
                              f"query string or confirm search manually "
                              f"at ENA: {url}")

        genome_sizes = get_genome_sizes()

        sanitized_dict = {}
        # Iterate over rows of dataframe to sanitize query entries:
        for index, entry in df.iterrows():

            # FTP Links
            ftp_links = entry["fastq_ftp"].strip(";").split(";")
            ftp_sizes = entry["fastq_bytes"].strip(";").split(";")

            if entry["library_layout"] == "PAIRED":

                # If there are fewer than two links, ignore accession:
                # this need to be fixed see project accession PRJEB12470
                # for ST59

                # Need to split FASTQ into forward and reverse after download:
                if len(ftp_links) < 2:
                    continue
                else:
                    # Get last two links and file sizes,
                    # which should be forward + reverse,
                    # check that they are conforming to
                    # pattern _1.fastq.gz and _2.fastq.gz:
                    ftp_links = ftp_links[-2:]
                    if "_1.fastq.gz" not in ftp_links[0] and \
                            "_2.fastq.gz" not in ftp_links[1]:
                        continue
                    ftp_1, ftp_2 = ftp_links
                    ftp_sizes = ftp_sizes[-2:]

            elif entry["library_layout"] == "SINGLE":

                # If there are fewer than one or more
                # than two links, ignore accession:
                if len(ftp_links) != 1:
                    continue
                else:
                    ftp_1 = ftp_links[0]
                    ftp_2 = None
            else:
                raise ValueError("Layout must be either SINGLE or PAIRED")

            # Convert to MB:
            try:
                size = sum([int(byte)/1024/1024 for byte in ftp_sizes])
            except ValueError:
                size = None

            try:
                reads = int(entry["read_count"])
            except ValueError:
                reads = None

            try:
                bases = int(entry["base_count"])
            except ValueError:
                bases = None

            if bases:
                try:
                    coverage = bases/(float(
                        genome_sizes.loc[entry["tax_id"], "Size"]
                    )*1000000)
                except KeyError or ValueError or ZeroDivisionError:
                    stamp(
                        "TaxID", entry["tax_id"],
                        " for determination of genome size could not be found."
                    )
                    coverage = None
            else:
                coverage = None

            entry_dict = {
                "id": str(uuid.uuid4()),
                "ftp_1": ftp_1,
                "ftp_2": ftp_2,
                "size": size,
                "reads": reads,
                "bases": bases,
                "coverage": coverage,
                "layout": entry["library_layout"],
                "platform": entry["instrument_platform"],
                "model": entry["instrument_model"],
                "source": entry["library_source"],
                "strategy": entry["library_strategy"],
                "tax_id": entry["tax_id"],
                "sample": entry["sample_accession"],
                "study": entry["study_accession"]
            }

            sanitized_dict[entry["run_accession"]] = entry_dict

        try:
            df = pandas.DataFrame(sanitized_dict).T
        except EmptyDataError:
            raise SurveyError(f"No results were returned for query: {url}")

        return df

    @staticmethod
    def _construct_species_query(
            species, platform, source, layout, strategy
    ) -> str:
        """ Construct a specific ENA query string to obtain paired-end
        Illumina reads or Nanopore reads from the DataWarehouse. """

        return f'tax_name("{species}") AND instrument_platform={platform}' \
            f' AND library_source={source} AND library_layout={layout}' \
            f' AND library_strategy={strategy}'

    @staticmethod
    def _construct_study_query(study):
        """Construct a specific ENA query
        string to obtain reads from a Study Accession"""

        return f'"study_accession={study}"'

    @staticmethod
    def _construct_sample_query(sample):

        return f'"sample_accession={sample}"'
