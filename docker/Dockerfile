# Couldn't get github actions to cache the docker builds so this starts from a sgpy_base image
FROM chasemc2/sgpy-base:1.0.0

ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python will not be found)

RUN /opt/conda/bin/pip install -i https://test.pypi.org/simple/ socialgene[full]==2.5.2 --extra-index-url https://pypi.org/simple

ENV PATH="/opt/conda/bin/neo4j/bin:$PATH:/opt/conda/condabin:/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
ENV NEO4J_HOME='/opt/conda/bin/neo4j'
ENV JAVA_HOME='/opt/conda'

WORKDIR $NEO4J_HOME

