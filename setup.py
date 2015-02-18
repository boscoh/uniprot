#!/usr/bin/env python
from setuptools import setup

description = "Docs at http://github.com/boscoh/uniprot"

setup(
    name='uniprot',
    version='1.1',
    author='Bosco Ho',
    author_email='boscoh@gmail.com',
    url='http://github.com/boscoh/uniprot',
    description='retrieve protein sequence identifiers and metadata from http://uniprot.org',
    long_description=description,
    license='BSD',
    install_requires=['requests'],
    py_modules=['uniprot'],
    scripts=['seqidtype']
)