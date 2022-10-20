include common_parameters.env
export

$(if $(wildcard ../sgpy),,$(error Why aren't you in the sgpy directory?))

SHELL = /bin/sh
CURRENT_UID := $(shell id -u)
CURRENT_GID := $(shell id -g)
export CURRENT_UID
export CURRENT_GID


## help	:	Print help info
help : Makefile
	@sed -n 's/^##//p' $<

## clean	:	Clean temporary files, INCLUDES DELETING THE NEXTFLOW TEMPORY WORK DIRECTORY
clean:
	@echo "Cleaning up temporary Python files..."
	find . -type f -name "*.py[co]" -exec rm -r {} +
	find . -type d -name "__pycache__" -exec rm -r {} +
	find . -type d -name "htmlcov" -exec rm -r {} +
	find . -type d -name "Autometa.egg-info" -exec rm -r {} +
	find . -type d -name "dist" -exec rm -r {} +
	find . -type d -name "build" -exec rm -r {} +
	find . -type d -name ".eggs" -exec rm -r {} +
	find . -type d -name "*.egg-info" -exec rm -r {} +
	find . -type d -name ".pytest_cache" -exec rm -r {} +
	find . -type d -name ".tox" -exec rm -r {} +
	find . -name ".coverage" -exec rm -r {} +
	find . -name "coverage.xml" -exec rm -r {} +
	find . -name ".nextflow.log.*" -exec rm -r {} +
	find . -type d -name ".nextflow" -exec rm -r {} +
	find . -name ".nextflow.log" -exec rm {} +
	find . -type d -name "work" -exec rm -r {} +
	find . -type d -name ".pytest_cache" -exec rm -r {} +
	find . -type d -name "__pycache__" -exec rm -r {} +
	find . -type d -name "site" -exec rm -r {} +
	find . -type d -name ".tox" -exec rm -r {} +


## socialgene_image	:	Make socialgene python docker image
build_docker_image: clean
	docker build -f Dockerfile --tag socialgene:latest .

build_docker_image2: clean
	docker build -f Dockerfileantismash --tag antismashsg:latest .

## conda :	Use conda to create the socialgene conda enviornment
conda:
	conda env create --file environment.yml

## mamba :	Use mamba to create the socialgene conda enviornment
mamba:
	mamba env create --file environment.yml

## install :	Install the socialgene python package 
install: 
	pip3 install -e .

pipinstall_test_extras:
	python3 -m pip install --upgrade pip setuptools wheel numpy pytest-cov

## pytest	:	Run Python pacakge unit tests
pytest: clean install pipinstall_test_extras 
	pytest tests -v --ignore=socialgene/entrypoints/export_protein_loci_assembly_tables.py 	 --cov=./socialgene --cov-report=xml:./coverage.xml --cov-report html
	xdg-open htmlcov/index.html
	
## pytestnf :	Run Nextflow pytest tests (first runs clean, install python and  nextflow test run)
pytestnf: clean install pipinstall_test_extras testnf
	python3 -m pip install --upgrade pip setuptools wheel numpy pytest-cov
	coverage run --source=./socialgene --module pytest ./socialgene/tests/nextflow --neo4j_outdir $(neo4j_outdir) 

run_ci: clean install pytest
	@echo "TESTING WITH FLAKE8"
	flake8 . --ignore E203,E501,W503 --count --show-source
	@echo "TESTING WITH BLACK"
	black --check . --extend-exclude temp_scripts

## run_ci :	Uploda to PyPi
run_cd: run_ci clean
	python3 -m pip install --user --upgrade setuptools wheel twine
	python3 -m build ./socialgene
	python3 -m twine upload --repository testpypi ./socialgene/dist/*

# https://stackoverflow.com/a/6273809/1826109
%:
	@:
###########################################################################################################
###########################################################################################################
