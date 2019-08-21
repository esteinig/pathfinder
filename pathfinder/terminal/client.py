import click

from .download import download
from .collect_survey import collect_survey
from .filter_survey import filter_survey


VERSION = '0.2'

@click.group()
@click.version_option(version=VERSION)
def terminal_client():
    pass


terminal_client.add_command(download)
terminal_client.add_command(collect_survey)
terminal_client.add_command(filter_survey)