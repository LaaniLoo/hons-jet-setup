# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import subprocess, os

# -- Project information -----------------------------------------------------

project = 'Jet Setup'
copyright = '2020, Patrick Yates-Jones'
author = 'Patrick Yates-Jones'


# -- General configuration ---------------------------------------------------

# Sphinx Multiversion
smv_srcdir = os.environ.get('SPHINX_MULTIVERSION_SOURCEDIR', None)
smv_tag_whitelist = r'^.*$'
smv_remote_whitelist = r'^(origin)$'
smv_released_pattern = r'^tags/.*$'
# smv_prefer_remote_refs = True

# Doxygen directory setup
if smv_srcdir is None:
    doxygen_dir = './'
    proj_dir = '../xml'
else:
    doxygen_dir = f'{smv_srcdir}/../'
    proj_dir = f'{smv_srcdir}/../xml'


# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "breathe",
    "sphinx_multiversion",
    # "exhale",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme_path = ['.']
html_theme = '_theme'
html_theme_options = {
    'collapse_navigation': False,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

breathe_projects = {
    "JetSetup":proj_dir,
}

breathe_default_project = "JetSetup"

breathe_domain_by_extension = {
    'h': 'c',
    'c': 'c',
}

exhale_args = {
    "containmentFolder": "./api",
    "rootFileName": "library_root.rst",
    "rootFileTitle": "Library API",
    "doxygenStripFromPath": "../../src/",
    "createTreeView": True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin": "INPUT = ../../src",
    "verboseBuild": True,
}

primary_domain = 'c'
highlight_language = 'c'

def run_doxygen(folder):
    """Run the doxygen make command in the designated folder"""
    try:
        retcode = subprocess.run('doxygen', shell=True, cwd=folder).returncode
        if retcode < 0:
            sys.stderr.write(f"doxygen terminated by signal {-retcode}")
    except OSError as e:
        sys.stderr.write(f"doxygen execution failed: {e}")

def generate_doxygen_xml(app):
    """Run the doxygen make commands"""

    run_doxygen(doxygen_dir)

def setup(app):
    # Add hook for building doxygen XML when needed
    app.connect("builder-inited", generate_doxygen_xml)
