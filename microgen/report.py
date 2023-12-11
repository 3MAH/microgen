"""Module managing errors."""
from collections import namedtuple
import re
import subprocess
import sys
import warnings

import scooby

class Report(scooby.Report):
    """Generate a microgen software environment report.

    Parameters
    ----------
    additional : sequence[types.ModuleType], sequence[str]
        List of packages or package names to add to output information.

    ncol : int, default: 3
        Number of package-columns in html table; only has effect if
        ``mode='HTML'`` or ``mode='html'``.

    text_width : int, default: 80
        The text width for non-HTML display modes.

    sort : bool, default: False
        Alphabetically sort the packages.

    Examples
    --------
    >>> import microgen
    >>> microgen.Report()
    --------------------------------------------------------------------------------
      Date: Thu Nov 23 16:33:48 2023 CET

                    OS : Darwin
                CPU(s) : 12
              Machine : arm64
          Architecture : 64bit
          Environment : IPython

      Python 3.11.6 | packaged by conda-forge | (main, Oct  3 2023, 10:37:07)
      [Clang 15.0.7 ]

              pyvista : 0.42.2
                  vtk : 9.2.5
                numpy : 1.26.0
                scooby : 0.7.4
                  gmsh : 4.11.1
              cadquery : 2.3.1
            matplotlib : 3.8.0
            pytest-cov : 4.1.0
                meshio : 5.3.4
    --------------------------------------------------------------------------------
    """

    def __init__(self, additional=None, ncol=3, text_width=80, sort=False):
        """Generate a :class:`scooby.Report` instance."""

        # Mandatory packages
        core = ['pyvista', 'vtk', 'numpy', 'scooby', 'gmsh', 'cadquery']

        # Optional packages.
        optional = [
            'matplotlib',
            'pytest-cov',
            'meshio'
        ]

        scooby.Report.__init__(
            self,
            additional=additional,
            core=core,
            optional=optional,
            ncol=ncol,
            text_width=text_width,
            sort=sort,
        )
