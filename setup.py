from setuptools import setup, find_packages

setup(
    name="pathfinder",
    url="https://github.com/pf-core/pathfinder",
    author="Eike J. Steinig",
    author_email="eikejoachim.steinig@my.jcu.edu.au",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "tqdm",
        "pandas",
        "colorama",
        "delegator.py",
        "google-cloud-storage",
        "click",
        "pytest",
        "pre-commit",
    ],
    entry_points="""
        [console_scripts]
        pf=pathfinder.terminal.client:terminal_client
        pathfinder=pathfinder.terminal.client:terminal_client
    """,
    version="0.2",
    license="MIT",
    description="Reproducible research and population genomic analysis "
                "pipelines for nanopore and Illumina data; genome surveillance "
                "of antimicrobial resistance and bacterial pathogens",
)
