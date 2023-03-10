FROM mambaorg/micromamba:1.1.0


COPY . .

RUN micromamba install -y -n base -f environment.yaml && \
    micromamba clean --all --yes

### Add Neo4j, remove then re-add the needed directories/files to avoid any potential permission issues
RUN /opt/conda/bin/wget -q https://neo4j.com/artifact.php?name=neo4j-community-5.1.0-unix.tar.gz -O neo4j-community-5.1.0-unix.tar.gz \
    && tar -xf neo4j-community-5.1.0-unix.tar.gz  \
    && rm neo4j-community-5.1.0-unix.tar.gz \
    && chmod -R 777 neo4j-community-5.1.0  \
    && mkdir -p /opt/conda/bin/neo4j \
    && mv neo4j-community-5.1.0  /opt/conda/bin/neo4j \
    && rm -rf /opt/conda/bin/neo4j/data /opt/conda/bin/neo4j/logs /opt/conda/bin/neo4j/import \
    && mkdir -p /opt/conda/bin/neo4j/data /opt/conda/bin/neo4j/logs /opt/conda/bin/neo4j/import \
    && touch /opt/conda/bin/neo4j/import.report \
    && chmod -R 777 /opt/conda/bin/neo4j/data /opt/conda/bin/neo4j/logs /opt/conda/bin/neo4j/import /opt/conda/bin/neo4j/import.report


USER $MAMBA_USER
ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python will not be found)

RUN /opt/conda/bin/pip3 install -e .

ENV PATH="/opt/conda/envs/bin:$PATH:/opt/conda/condabin:/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/conda/bin/neo4j/neo4j-community-5.1.0/bin"

WORKDIR /opt/conda/bin/neo4j/neo4j-community-5.1.0
