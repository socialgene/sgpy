Note: This is all pre-alpha stuff (i.e. being worked on extensively, there will be breaking changes, the repo may be burnt down and rebuilt at any time). Extensive documentation will be made available at a later date when this is ready for general use.

![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/socialgene/sgpy)
[![Linting](https://github.com/socialgene/sgpy/actions/workflows/linters.yml/badge.svg)](https://github.com/socialgene/sgpy/actions/workflows/linters.yml)
[![Continuous Integration](https://github.com/socialgene/sgpy/actions/workflows/pr_ci.yml/badge.svg)](https://github.com/socialgene/sgpy/actions/workflows/pr_ci.yml)
[![Continuous Deployment](https://github.com/socialgene/sgpy/actions/workflows/pypi_autodeploy_python.yml/badge.svg)](https://github.com/socialgene/sgpy/actions/workflows/pypi_autodeploy_python.yml)

Documentation can be found here: <https://socialgene.github.io>

<!---
To create the UML diagram of the library:
```bash
pyreverse -o png -p sgpy socialgene
```
--->

![classes](classes_sgpy.png)
![packages](packages_sgpy.png)





## Design

The code is organized under a number of submodules/directories:

- base: core functions of the library
- cli: all command line interface code
- clustermap: used to convert a socialgene object to clustermap json
- findmybgc
- hashing
- hmm: code for working with HMMER
- neo4j: code for working with SocialGene Neo4j databases
- parsers: external file parsers (e.g. genbank, fasta, HMMER results, etc)
- scoring: functions for measuring protein similarity
- taxonomy
- utils

## Installation with pip

<https://pypi.org/project/socialgene>

```bash
pip install socialgene
```

## Create conda environment and install python package inside

```bash
git clone https://github.com/socialgene/sgpy.git
cd sgpy
make create_conda
```

## Build Python package from source

```bash
git clone https://github.com/socialgene/sgpy.git
cd sgpy
make install_python
```

## Build local Docker image

```bash
git clone https://github.com/socialgene/sgpy.git
cd sgpy
make build_docker_image
```


## Run pytest tests

```bash
git clone https://github.com/socialgene/sgpy.git
cd sgpy
make create_conda
make pytest
```

## Run all tests

```bash
git clone https://github.com/socialgene/sgpy.git
cd sgpy
make create_conda
make run_ci
```

## User-facing classes

### `SocialGene()`

This is the main class that most other user-facing classes should/do inherit from

### `FindMyBGC()`

### `SingleProteinSearch()`

#### Common example use cases

Starting with a single input protein and

- [want to compare it against all other proteins in the Neo4j database](jupyter/single_protein_search.ipynb)

Starting with a set of proteins (BGC) and

- [want to compare against all other proteins in the Neo4j database](jupyter/findmybgc.ipynb)

## Other

Most of the the classes that describe the structure of `SocialGene()` (e.g. proteins, domains, loci) live in `socialgene/src/socialgene/classes/molbio.py`
