FROM mambaorg/micromamba:1.1.0

COPY environment.yaml /tmp/env.yaml

RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

### Add Neo4j, remove then re-add the needed directories/files to avoid any potential permission issues
# https://neo4j.com/docs/operations-manual/current/configuration/file-locations/
RUN /opt/conda/bin/wget -q https://neo4j.com/artifact.php?name=neo4j-community-5.7.0-unix.tar.gz -O neo4j.tar.gz \
    && tar -xf neo4j.tar.gz \
    && rm neo4j.tar.gz \
    && mv neo4j* neo4j \
    && chmod -R 777 neo4j  \
    && mkdir -p /opt/conda/bin/neo4j \
    && mv neo4j  /opt/conda/bin/ \
    && rm -rf /opt/conda/bin/neo4j/data /opt/conda/bin/neo4j/logs /opt/conda/bin/neo4j/import \
    && mkdir -p /opt/conda/bin/neo4j/data /opt/conda/bin/neo4j/logs /opt/conda/bin/neo4j/import \
    && touch /opt/conda/bin/neo4j/import.report \
    && chmod -R 777 /opt/conda/bin/neo4j/data /opt/conda/bin/neo4j/logs /opt/conda/bin/neo4j/import /opt/conda/bin/neo4j/import.report

ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python will not be found)

ENV PATH="/opt/conda/bin/neo4j/bin:$PATH:/opt/conda/condabin:/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
ENV NEO4J_HOME='/opt/conda/bin/neo4j'

ENV JAVA_HOME='/opt/conda'

WORKDIR $NEO4J_HOME

