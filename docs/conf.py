# Configuration file for the Sphinx documentation builder.
#
# Theme and gallery set-up follow the fedoo and simcoon documentations, so
# that the three 3MAH packages share the same look.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import datetime
import os
import sys

import pyvista
from pyvista.plotting.utilities.sphinx_gallery import DynamicScraper
from sphinx_gallery.sorting import ExplicitOrder, FileNameSortKey

sys.path.insert(0, os.path.abspath(".."))

import microgen

# -- Project information -----------------------------------------------------

project = "microgen"
author = "3MAH development team"
copyright = f"2022-{datetime.date.today().year}, {author}"
release = microgen.__version__
version = release

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "jupyter_sphinx",
    "myst_parser",
    "sphinx_design",
    "sphinx_gallery.gen_gallery",
    "pyvista.ext.plot_directive",
    "pyvista.ext.viewer_directive",
]

templates_path = ["_templates"]
master_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_logo = "_static/microgen.png"
html_title = f"microgen {release}"
html_theme_options = {
    "canonical_url": "https://3mah.github.io/microgen-docs/",
    "logo_only": True,
    "style_external_links": True,
    "style_nav_header_background": "#24445C",
    "collapse_navigation": False,
    "navigation_depth": 3,
}
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_show_sphinx = False

# -- PyVista and Sphinx-Gallery ----------------------------------------------

pyvista.OFF_SCREEN = True
pyvista.BUILDING_GALLERY = True
pyvista.set_plot_theme("document")
pyvista.global_theme.window_size = [1024, 768]

# Interactive scenes heavier than this are replaced by the static image, to
# keep the published documentation small (dense TPMS meshes reach 100 MB).
MAX_INTERACTIVE_SCENE_BYTES = 5_000_000


class SizeCappedScraper(DynamicScraper):
    """PyVista scraper with interactive scenes only for light meshes."""

    def __call__(self, block, block_vars, gallery_conf):
        for plotter in pyvista.plotting.plotter._ALL_PLOTTERS.values():
            vtksz = getattr(plotter, "last_vtksz", None)
            if vtksz is not None and len(vtksz) > MAX_INTERACTIVE_SCENE_BYTES:
                plotter.last_vtksz = None
        return super().__call__(block, block_vars, gallery_conf)


# Examples that are not run at each build: too slow (more than 50 s) or
# depending on the mmg executables, which are not installable with pip. Their
# page shows the code and a precomputed image.
_not_executed = [
    "3Doperations/voronoi_gyroid.py",
    "Mesh/mmg3d.py",
    "Mesh/mmg_voronoi.py",
    "Mesh/remesh.py",
    "Mesh/periodic_gyroid_remesh.py",
]

sphinx_gallery_conf = {
    "examples_dirs": ["../examples/"],
    "gallery_dirs": ["examples"],
    "filename_pattern": r"^(?!.*/(?:"
    + "|".join(path.replace(".", r"\.") for path in _not_executed)
    + r")$).*\.py$",
    "ignore_pattern": r"run_examples\.py",
    "subsection_order": ExplicitOrder(
        [
            f"../examples/{category}"
            for category in (
                "BasicShapes",
                "TPMS",
                "Spinodoids",
                "Lattices",
                "ImplicitShapes",
                "RVE",
                "3Doperations",
                "Mesh",
                "Fibers",
            )
        ]
    ),
    "within_subsection_order": FileNameSortKey,
    # One toctree per section of the examples index, so the side bar lists
    # each example under its section (as in the simcoon documentation).
    "nested_sections": False,
    "image_scrapers": (SizeCappedScraper(), "matplotlib"),
    "download_all_examples": False,
    "remove_config_comments": True,
}
