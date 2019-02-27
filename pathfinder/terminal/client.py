import click

from .collect import collect
from .get import get
from .run import run
from .store import store
from .survey import survey
from .setup import setup
from .config import config
from .meta import meta
from .app import app
from .ref import ref

VERSION = '0.2'


@click.group()
@click.version_option(version=VERSION)
def terminal_client():
    pass


terminal_client.add_command(collect)
terminal_client.add_command(get)
terminal_client.add_command(run)
terminal_client.add_command(config)
terminal_client.add_command(store)
terminal_client.add_command(survey)
terminal_client.add_command(setup)
terminal_client.add_command(meta)
terminal_client.add_command(app)
terminal_client.add_command(ref)