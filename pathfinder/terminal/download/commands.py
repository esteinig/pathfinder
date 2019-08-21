import click
import pandas

from pathfinder.survey import Survey
from pathfinder.survey import MiniAspera

from pathlib import Path


@click.command()
@click.option(
    '--outdir', '-o', type=str, default="pf-download",
    help='Output directory for read files'
)
@click.option(
    '--batch', '-b', type=int, default=0,
    help='Batch large search results into their own output directories'
)
@click.option(
    '--file', '-f', type=str, default=0,
    help='CSV file with column: accession or project, to download.'
)
@click.option(
    '--accession', '-a', type=str, default=None,
    help='Accession or comma-separated string of accessions'
)
@click.option(
    '--project', '-p', type=str, default=None,
    help='Project accession or comma-separated string of project accessions'
)
@click.option(
    '--species', '-s', type=str, default=None,
    help='Scientific species name or TaxID'
)
@click.option(
    '--query', '-q', type=str, default=None,
    help='Custom search query string for ENA warehouse or '
         'path to CSV file from previous query.'
)
@click.option(
    '--filter', '-f', type=str, default=None,
    help='Custom query filter on the return fields of the ENA query '
         'for example: "coverage < 700 & coverage > 50" '
)
@click.option(
    '--scheme', type=str, default=None,
    help='Sequence read scheme, either: Illumina (WGS, PE, GENOMIC), '
    'Nanopore (WGS, SINGLE, GENOMIC) or None (default) '
    'for filtering ENA query results by type templates'
)
@click.option(
    '--ftp', '-f', is_flag=True,
    help='Force download from FTP instead of Aspera (slow)'
)
@click.option(
    '--limit', '-l', type=int, default=0,
    help='Limit download to first --limit query results.'
)
@click.option(
    '--submitted', is_flag=True,
    help='Use default FASTQ files from ENA, switch on to use '
         'read files uploaded by submitter.'
)
def download(
    outdir,
    batch,
    file,
    accession,
    project,
    species,
    query,
    filter,
    scheme,
    ftp,
    limit,
    submitted
):
    """ Download sequence read data from ENA """

    if file:
        df = pandas.read_csv(file)
        df.columns = [c.lower() for c in df.columns]
        if 'accession' in df:
            accession = df['accession'].tolist()
            project, species, query = (None,) * 3
        elif 'project' in df:
            project = df['project']
            accession, species, query = (None,) * 3
        else:
            print('Could not find columns: accession or project in CSV.')

    Path(outdir).mkdir(exist_ok=True, parents=True)

    survey = Survey(outdir=outdir)

    if query is not None:
        if Path(query).exists():
            query_csv = Path(query)
            survey.query_from_csv(query_csv)
    else:
        survey.query_ena(
            sample=accession,
            study=project,
            species=species,
            term=query,
            scheme=scheme,
            submitted_fastq=submitted
        )

        print(survey.query)

        query_csv = Path(f"{outdir}/query.csv")
        survey.query_to_csv(query_csv)

    if filter is not None:
        survey.filter_query(filter)

    if batch > 0:
        batches = survey.batch(batch_size=batch)
        batches = survey.batch_output(
            batches, outdir=Path(outdir)
        )
    else:
        batches = [(
            Path(outdir), query_csv
        )]

    for batch_path, batch_csv in batches:
        ascp = MiniAspera()
        ascp.download_batch(
            file=batch_csv,
            outdir=batch_path,
            limit_download=limit,
            ftp=ftp
        )