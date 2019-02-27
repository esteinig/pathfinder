#!/usr/bin/env python

""" Class for pipeline management in Pathfinder.

This module contains the PipelineManager class, which
exposes methods to obtain current information about the
available pipelines from the Github repository and enables
users to pull / setup pipelines in local storage.

"""

from pathlib import Path
from github import Github

# Subclass Pipeline in specific objects


class Pipeline:
    """ Pipelines base object to store and test pipeline architectures """
    def __init__(self, name, version, license, prefix='pf-'):
        """ Object for pipeline repository tree structure -
        align with nf-core for transferability when we make
        the pipelines open source:

            |____assets
            |____bin
            |____conf
            |__test.config
            |____docs
            |__.gitignore
            |__.travis.yml
            |__CHANGELOG.md
            |__Dockerfile
            |__LICENSE.md
            |__README.md
            |__Singularity
            |__environment.yml
            |__main.nf
            |__nextflow.config

        """

        self.name = name
        self.version = version
        self.license = license

        self.readme = None
        self.changelog = None

        self.gitignore = None

        self.docker = None
        self.conda = None
        self.singularity = None

        self.prefix = prefix

    # Test accessors

    def has_name(self):
        """ Does the pipeline have a name and does it start with prefix? """
        return True if self.name and self.name.startswith(self.prefix) else False

    def has_version(self):
        """ Does the pipeline have a version tag ? """
        return True if self.version else False

    def has_license(self):
        """ Does the pipeline have a license name? """
        return True if self.license else False

    def has_readme(self):
        """ Does the pipeline have a MREADME.md ? """
        return True if self.readme else False

    def has_gitignore(self):
        """ Does the pipeline have a .gitignore file? """
        return True if self.gitignore else False

    def has_conda(self):
        """ Does the pipeline have an environment.yaml file? """
        return True if self.gitignore else False

    def has_docker(self):
        """ Does the pipeline have a Dockerfile? """
        return True if self.docker else False


class PipelineManager:

    def __init__(self, user='pf-core', prefix='pf-'):

        self.github = Github()
        self.repositories = self._get_repos(user=user, prefix=prefix)

        self.pipelines = dict()

        self.nextflow = Path().home() / '.nextflow'

    def get_info(self):
        """Extract information from pipeline repositories"""
        for repo in self.repositories:
            tags = list(repo.get_tags())
            license_file = repo.get_license()

            try:
                license_name = license_file.license.name
            except AttributeError:
                license_name = ''
            try:
                pipeline_name = repo.name
            except AttributeError:
                pipeline_name = ''

            self.pipelines[pipeline_name] = {
                'tag': tags[0].name if tags else '',
                'license': license_name,
                'repo': repo
            }

    def print_info(self):
        """ Print a table with information on available pipelines, versions and licenses."""
        print(f"{'Pipeline':<25}{'Version':<20}{'License':<25}")
        print(11*'-'+14*' '+11*'-'+9*' '+11*'-'+14*' ')
        for pipeline in sorted(self.pipelines):
            print(f"{pipeline:<25}{self.pipelines[pipeline]['tag']:<20}{self.pipelines[pipeline]['license']:<25}")

    def _get_repos(self, user, prefix):
        """Get repositories from connection"""
        return [repo for repo in self.github.search_repositories(f'user:{user}') if prefix in repo.name]


pp = PipelineManager()
pp.get_info()
pp.print_info()

# Rules: fstring should be "" so that '' can be used inside
# single line return functions comments without newline, otherwise newline
# pipeline names < 25 charachters, versio ntags < 15 characters, license names <