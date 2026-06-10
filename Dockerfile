FROM mambaorg/micromamba:1.5.8

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# install the python package itself
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/fragmentomics
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install --no-deps /opt/fragmentomics

ENV PATH=/opt/conda/bin:$PATH
