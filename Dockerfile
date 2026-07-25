FROM continuumio/miniconda3:latest

LABEL maintainer="Brown Beckley <brownbeckley94@gmail.com>"
LABEL description="PseudoScope - Complete P. aeruginosa genomic analysis pipeline"

RUN apt-get update && apt-get install -y \
    procps \
    jq \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/pseudoscope

COPY . /opt/pseudoscope/

ENV CONDA_SHARDED_REPODATA=false

RUN conda env create -f environment.yml && \
    conda clean -afy

SHELL ["conda", "run", "-n", "pseudoscope", "/bin/bash", "-c"]

RUN abricate --setupdb

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "pseudoscope", "pseudoscope"]
CMD ["-h"]
