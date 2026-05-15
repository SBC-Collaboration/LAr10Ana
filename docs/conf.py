# Configuration file for the Sphinx documentation builder.
import os, sys

# -- Project information

project = 'SBC LAr10Ana'
copyright = '2024-2026 SBC Collaboration'
author = 'SBC Collaboration'

release = '0.0.3'
version = '0.0.3'

# -- General configuration

extensions = [
    "myst_parser",
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

root_doc = 'index'

# Add icon
html_static_path = ['assets']
html_favicon = 'assets/sbc_icon.png'

# Set each section is a chapter in pdf file
latex_toplevel_sectioning = 'chapter'

# Remove blank pages in pdf file
latex_elements = {
  'extraclassoptions': 'openany,oneside'
}

sys.path.insert(0, os.path.abspath(".."))

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'
