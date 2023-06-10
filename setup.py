#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from glob import glob
from os.path import basename, splitext

############
#  https://packaging.python.org/en/latest/guides/making-a-pypi-friendly-readme/
from pathlib import Path

from setuptools import find_packages, setup

this_directory = Path(__file__).parent
long_description = Path(this_directory / "README.md").read_text()
############

setup(
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
    python_requires=">=3.11",
)
