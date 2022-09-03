FROM python:3.10-slim-bullseye

RUN apt-get update && apt-get install --no-install-recommends -y \
    # dependencies for building Python packages
    gcc \
    build-essential \
    # procps is required to be in docker images used by nextflow
    procps \
    # other things w might want within the nextflow pipeline
    wget \
    rsync \
    hmmer \
    pigz

RUN /usr/local/bin/python -m pip install --upgrade pip

RUN pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple socialgene

