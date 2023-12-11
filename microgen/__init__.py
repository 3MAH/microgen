import importlib.metadata

from .external import *
from .mesh import *
from .operations import *
from .periodic import *
from .phase import *
from .report import *
from .rve import *
from .shape import *
from .test_utils import *

__version__ = importlib.metadata.version(__package__ or __name__)
