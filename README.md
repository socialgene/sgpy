Note: All releases may intorduce breaking changes until the release of v1.0.0

![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/socialgene/sgpy)
[![codecov](https://codecov.io/gh/socialgene/sgpy/branch/main/graph/badge.svg?token=8f8GCc4J3G)](https://codecov.io/gh/socialgene/sgpy)
[![Linting](https://github.com/socialgene/sgpy/actions/workflows/linters.yml/badge.svg)](https://github.com/socialgene/sgpy/actions/workflows/linters.yml)
[![Continuous Integration](https://github.com/socialgene/sgpy/actions/workflows/pr_ci.yml/badge.svg)](https://github.com/socialgene/sgpy/actions/workflows/pr_ci.yml)
[![Continuous Deployment](https://github.com/socialgene/sgpy/actions/workflows/pypi_autodeploy_python.yml/badge.svg)](https://github.com/socialgene/sgpy/actions/workflows/pypi_autodeploy_python.yml)

Both user and developer documentation can be found at: <https://socialgene.github.io>

![classes](classes_sgpy.png)
![packages](packages_sgpy.png)

<!---
To create the UML diagram of the library:
```bash
pyreverse -o png -p sgpy socialgene
```
--->

## Installation with pip
![https://pypi.org/project/socialgene](https://img.shields.io/pypi/dm/socialgene)

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

# Run all tests

```bash
git clone https://github.com/socialgene/sgpy.git
cd sgpy
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
