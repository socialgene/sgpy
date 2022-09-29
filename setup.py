#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from glob import glob
from os.path import basename
from os.path import splitext
from setuptools import find_packages
from setuptools import setup

############
#  https://packaging.python.org/en/latest/guides/making-a-pypi-friendly-readme/
from pathlib import Path

this_directory = Path(__file__).parent
long_description = Path(this_directory / "README.md").read_text()
############

setup(
    name="socialgene",
    version="1.0.1",
    license="MIT",
    description="Creating and interacting with graph databases of protein domains and their genome coordinates",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Chase M. Clark",
    author_email="chasingmicrobes@gmail.com",
    url="https://github.com/socialgene/sgpy",
    packages=find_packages(),
    package_data={
        "": [
            "*.env",
            "data/biosample_attributes",
            "neo4j/queries.cypher",
        ]
    },
    include_package_data=True,
    py_modules=[splitext(basename(path))[0] for path in glob("src/*.py")],
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: Implementation :: CPython",
        "Programming Language :: Python :: Implementation :: PyPy",
        "Topic :: Utilities",
    ],
    project_urls={
        "Documentation": "https://socialgene.github.io",
        "Changelog": "https://github.com/socialgene/sgpy/blob/main/CHANGELOG.md",
        "Issue Tracker": "https://github.com/socialgene/sgpy/issues",
    },
    keywords=[
        # eg: 'keyword1', 'keyword2', 'keyword3',
    ],
    python_requires=">=3.9",
    install_requires=[
        "rich>=10.12.0",
        "pandas>=1.3",
        "numpy>=1.21",
        "neo4j>=4.3",
        "biopython>=1.7",
        "textdistance>=4.2.1",
    ],
    entry_points={
        "console_scripts": [
            "socialgene = socialgene.cli.__main__:main",
            "socialgene_process_domtblout = socialgene.cli.process_domtblout:main",
            "socialgene_clean_hmm = socialgene.cli.clean_hmms:main",
            "socialgene_process_genbank= socialgene.cli.export_protein_loci_assembly_tables:main",
            "socialgene_ncbi_taxonomy= socialgene.taxonomy.ncbi_taxonomy:main",
            "socialgene_export_parameters = socialgene.cli.parameter_export:main",
            "socialgene_export_neo4j_headers = socialgene.cli.export_neo4j_header_files:main",
            "socialgene_create_neo4j_db = socialgene.cli.create_neo4j_db:main",
            "socialgene_version = socialgene.utils.version:main",
            "socialgene_hmm_tsv_parser= socialgene.cli.socialgene_hmm_tsv_parser:main",
            "socialgene_prothash_sqlite= socialgene.cli.protein_sqlite:main",
        ]
    },
)
