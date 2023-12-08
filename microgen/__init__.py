import importlib.metadata

from .external import *
from .mesh import *
from .operations import *
from .periodic import *
from .phase import *
from .rve import *
from .shape import *
from .report import *

__version__ = importlib.metadata.version(__package__ or __name__)
