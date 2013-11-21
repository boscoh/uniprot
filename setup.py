#!/usr/bin/env python
from setuptools import setup


description = \
"""uniprot.py is a collection of calls to retrieve protein sequence identifier mappings and metadata from the uniprot.org website.

Docs at http://github.com/boscoh/uniprot.
"""

setup(
    name='uniprot',
    version='1.0a',
    author='Bosco Ho',
    author_email='boscoh@gmail.com',
    url='http://github.com/boscoh/uniprot',
    description='retrieve protein sequence identifiers and metadata from uniprot.org',
    long_description=description,
    license='GPLv3',
    install_requires=['requests'],
    py_modules=['uniprot',],
    scripts=['seqidtype.py']
)