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
