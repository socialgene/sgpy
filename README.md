Note: All releases may intorduce breaking changes until the release of v1.0.0

# Status

## Github

![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/socialgene/sgpy) [![codecov](https://codecov.io/gh/socialgene/sgpy/branch/main/graph/badge.svg?token=8f8GCc4J3G)](https://codecov.io/gh/socialgene/sgpy) [![Linting](https://github.com/socialgene/sgpy/actions/workflows/linters.yml/badge.svg)](https://github.com/socialgene/sgpy/actions/workflows/linters.yml) [![Continuous Integration](https://github.com/socialgene/sgpy/actions/workflows/pr_ci.yml/badge.svg)](https://github.com/socialgene/sgpy/actions/workflows/pr_ci.yml) [![Continuous Deployment](https://github.com/socialgene/sgpy/actions/workflows/pypi_autodeploy_python.yml/badge.svg)](https://github.com/socialgene/sgpy/actions/workflows/pypi_autodeploy_python.yml)
## PyPI

https://pypi.org/project/socialgene/

![PyPI](https://img.shields.io/pypi/v/socialgene) ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/socialgene) ![PyPI - Status](https://img.shields.io/pypi/status/socialgene) ![https://pypi.org/project/socialgene](https://img.shields.io/pypi/dm/socialgene)

# Contributing

Please see https://github.com/socialgene/sgpy/blob/main/CONTRIBUTING.md

# Documentation

Both user and developer documentation can be found at: <https://socialgene.github.io>

<!---
To create the UML diagram of the library:
```bash
pyreverse -o png -p sgpy socialgene
```
--->

## Installation with pip

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

## Run all tests

```bash
git clone https://github.com/socialgene/sgpy.git
cd sgpy
make run_ci
```


## Classes:
![classes](classes_sgpy.png)

## Modules
![modules](packages_sgpy.png)
