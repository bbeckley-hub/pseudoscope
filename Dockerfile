FROM continuumio/miniconda3:latest

LABEL maintainer="Brown Beckley <brownbeckley94@gmail.com>"
LABEL description="PseudoScope - Complete P. aeruginosa genomic analysis pipeline"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    procps \
    jq \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/pseudoscope

# Copy entire project
COPY . /opt/pseudoscope/

# Create Conda environment
RUN conda env create -f environment.yml && \
    conda clean -afy

# Make the environment the default for RUN commands
SHELL ["conda", "run", "-n", "pseudoscope", "/bin/bash", "-c"]

# Run abricate database setup (one-time)
RUN abricate --setupdb

# Run AMR database setup using the module's own script
RUN cd /opt/pseudoscope/pseudoscope/modules/amr_module && \
    python p_amrfinder.py --update-db

# Set entrypoint
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "pseudoscope", "pseudoscope"]
CMD ["-h"]
