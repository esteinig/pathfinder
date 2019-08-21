from pathlib import Path
import click

from pathfinder.pipelines.survey import SurveyResult


@click.command()
@click.option(
    '--data', '-d', type=Path, required=True,
    help='Directory containing data collected from pipeline.'
)
@click.option(
    '--outdir', '-o', type=Path, default="pf-survey-filtered",
    help='Output directory for filtered pipeline results.'
)
@click.option(
    '--species', type=str, default='Staphylococcus aureus',
    help='Species top hit filter for taxonomic identification: pf-survey'
)
@click.option(
    '--purity', type=float, default=0.8,
    help='Minimum threshold for taxonomic identification of species: pf-survey'
)
@click.option(
    '--contamination', type=float, default=0.02,
    help='Maximum allowed proportion of taxonomic hits '
         'belonging to other species'
)
@click.option(
    '--gene_coverage', type=float, default=0.8,
    help='Minimum coverage to retain typed genes with Abricate'
)
@click.option(
    '--gene_identity', type=float, default=0.9,
    help='Minimum identity to retain typed genes with Abricate'
)
@click.option(
    '--not_complete_lineage', is_flag=True,
    help='Include isolates with uncertainty in lineage assignment.'
)
@click.option(
    '--complete', '-c', is_flag=True,
    help='Output intersection of all genomes present'
         ' after filtering across analyses.'
)
def filter_survey(
        data,
        outdir,
        species,
        purity,
        contamination,
        gene_identity,
        gene_coverage,
        not_complete_lineage,
        complete
):
    """ Collect results from pipeline: pf-survey """

    sr = SurveyResult()

    sr.data.read(data)

    # Raw data can be written to file here if needed
    sr.filter(
        species=species,
        purity=purity,
        contamination=contamination,
        gene_coverage=gene_coverage,
        gene_identity=gene_identity,
        complete_lineage=not not_complete_lineage
    )

    if complete:
        sr.data.complete()

    sr.data.write(outdir)
