from dataclasses import dataclass
from collections import Callable

# Template for processes attribute in
# subclasses of `pathfinder.process.Processor`


@dataclass
class SurveyProcessSetting:
    """ Template for processes attribute in
    subclasses of `pathfinder.process.Processor` """

    files: list = None
    parse_func: Callable or None = None
    process_func: Callable or None = None
    remove: str = '.extension'
    aggregate: bool = True


@dataclass
class SurveyCleanerSetting:
    """ Convenience class for passing settings to
    SurveyResult and SurveyCleaner

    Default cleaner settings can be defined here
    with user access to change parameters in:
    `pathfinder.process.results.Result.clean()`
    which in turn inherits the methods from:
    `pathfinder.process.cleaners.Cleaner`
    """

    purity: float = 0.8
    species: int or str = None
    contamination: float = 0.02
    gene_identity: float = 0.9
    gene_coverage: float = 0.9
    complete_lineage: bool = True
