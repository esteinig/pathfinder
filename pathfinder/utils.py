import os
import sys
import time
import pathlib
import hashlib
import json
import pandas
import itertools
from pathlib import Path
from colorama import init, Style, Fore

# Colorama Initialization:

init()


def pretty_print(*args, color="yellow"):

    if color == "yellow":
        print(Fore.YELLOW + Style.BRIGHT + " ".join([str(arg) for arg in args]) + Style.RESET_ALL)
    elif color == "green":
        print(Fore.GREEN + Style.BRIGHT + " ".join([str(arg) for arg in args]) + Style.RESET_ALL)
    elif color == "red":
        print(Fore.RED + Style.BRIGHT + " ".join([str(arg) for arg in args]) + Style.RESET_ALL)


def stamp(*args, color="yellow"):

    if color == "yellow":
        print(str(time.strftime("[%H:%M:%S]")) + " " + Style.BRIGHT +
              Fore.YELLOW + " ".join([str(arg) for arg in args]) + Style.RESET_ALL)
    elif color == "green":
        print(str(time.strftime("[%H:%M:%S]")) + " " + Style.BRIGHT +
              Fore.GREEN + " ".join([str(arg) for arg in args]) + Style.RESET_ALL)
    elif color == "red":
        print(str(time.strftime("[%H:%M:%S]")) + " " + Style.BRIGHT +
              Fore.RED + " ".join([str(arg) for arg in args]) + Style.RESET_ALL)


def get_batches(it, size):

    """ Batch iterator from:
    https://stackoverflow.com/questions/28022223/
    how-to-iterate-over-a-dictionary-n-key-value-pairs-at-a-time
    """

    it = iter(it)
    while True:
        p = tuple(itertools.islice(it, size))
        if not p:
            break
        yield p


def get_simple_date(datetime_object):

    return datetime_object.strftime("%d-%m-%Y")


def get_genome_sizes():

    """
    Get median genome size for given species name from
    resources/prokaryote.sizes.txt (NCBI Prokaryot DB)
    by searching for TaxID.

    """

    genome_sizes = Path(__file__).parent / "resources" / "genome.sizes"
    genome_sizes = pandas.read_csv(genome_sizes, index_col=0)

    return genome_sizes


def get_package_path():

    """ Return directory path of file that initiates Python interpreter, i.e. main script for PathFinder """

    return Path(os.path.dirname(os.path.realpath(sys.argv[0])))


def get_project(path=os.getcwd()):

    """

    Walk up directory tree and attempt to find the first .pathfinder
    directory belonging to the current project. Return the project
    configuration directory, the configuration file (project.json)
    and database path (db/project.db).

    """

    path = pathlib.Path(path)

    for parent in [os.getcwd()] + list(path.parents):
        pathfinder = os.path.join(str(parent), ".pathfinder")
        print(pathfinder)
        if os.path.isdir(pathfinder):
            with open(os.path.join(pathfinder, "project.json")) as config_json:
                config = json.load(config_json)
            return pathfinder, config, os.path.join(pathfinder, "db", "project.db")

    return None, None, None


def md5(file):

    """ https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file """

    hash_md5 = hashlib.md5()
    with open(file, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def get_subdict(key, dictionary):

    for k, v in dictionary.items():
        if k == key:
            yield v
        elif isinstance(v, dict):
            for result in get_subdict(key, v):
                yield result
        elif isinstance(v, list):
            for d in v:
                if isinstance(d, dict):
                    for result in get_subdict(key, d):
                        yield result


# Pipeline Modules

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