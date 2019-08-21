from pathlib import Path
import click

from pathfinder.pipelines.survey import SurveyResult


@click.command()
@click.option(
    '--data', '-d', type=Path, required=True,
    help='Output directory containing result data from pipeline.'
)
@click.option(
    '--batch', '-b', is_flag=True,
    help='Process output of a batched run in directory of directories.'
)
@click.option(
    '--outdir', '-o', type=Path, default="pf-survey-result",
    help='Output directory for collection results.'
)
def collect_survey(
        data,
        outdir,
        batch,
):
    """ Collect and organise raw data from: pf-survey """

    outdir.mkdir(exist_ok=True, parents=True)

    if batch:
        dirs = [p / 'analysis' for p in list(
            Path(data).glob('batch_*')
        )]
    else:
        dirs = [Path(data)]

    survey_results = list()
    for path in dirs:
        print(f'Processing: {path}')
        sr = SurveyResult(path=path)
        sr.parse()

        # Something weird here, FASTA are added successively
        sr.data.get_file_paths(
            result_path=path, fasta_dir='skesa', fastq_dir='trimmomatic'
        )

        print('Fasta length', len(sr.data.fasta))

        if batch:
            sr.data.write(outdir / 'batch_data' / path.parent.stem)

        survey_results.append(sr)

    for result in survey_results:
        print(len(result.data.fasta))

    result = sum(survey_results)

    print(len(result.data.fasta))

    result.data.write(outdir)


