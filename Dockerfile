#syntax=docker/dockerfile:1.6

FROM mambaorg/micromamba:1.5.8-jammy

LABEL maintainer="Brown Beckley <brownbeckley94@gmail.com>"
LABEL description="PseudoScope - Complete P. aeruginosa genomic analysis pipeline"

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
        procps jq git ca-certificates wget unzip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/pseudoscope

COPY environment.lock /opt/pseudoscope/

RUN --mount=type=cache,target=/opt/conda/pkgs \
    --mount=type=cache,target=/root/.mamba/pkgs \
    micromamba create -y -n pseudoscope \
        --file /opt/pseudoscope/environment.lock \
        --no-rc --no-env

RUN --mount=type=cache,target=/opt/conda/pkgs \
    --mount=type=cache,target=/root/.mamba/pkgs \
    micromamba install -y -n pseudoscope \
        -c conda-forge -c bioconda \
        'prodigal>=2.6.3' \
        'diamond>=2.1.10' \
        'mlstdb>=1.3.0' \
        'fastani>=1.3.0'\
    --no-rc

COPY . /opt/pseudoscope/

RUN chmod -R a+rX /opt/pseudoscope/pseudoscope/modules/ && \
    find /opt/pseudoscope/pseudoscope/modules/ -type f \
        \( -name "*.py" -o -name "*.pl" -o -name "*.sh" \) \
        -exec chmod a+x {} \;

SHELL ["micromamba", "run", "-n", "pseudoscope", "/bin/bash", "-c"]

RUN pip install --no-deps --no-build-isolation -e /opt/pseudoscope

RUN abricate --setupdb

RUN if [ ! -d /opt/pseudoscope/pseudoscope/modules/amr_module/data/amrfinder_db ] || \
       [ -z "$(ls -A /opt/pseudoscope/pseudoscope/modules/amr_module/data/amrfinder_db 2>/dev/null)" ]; then \
        echo "[ERROR] AMR database missing from build context"; \
        echo "Expected: pseudoscope/modules/amr_module/data/amrfinder_db/20*/"; \
        exit 1; \
    fi

RUN MLST_DIR=/opt/pseudoscope/pseudoscope/modules/mlst_module/db/pubmlst/paeruginosa; \
    if [ ! -d "$MLST_DIR" ] || [ -z "$(ls -A "$MLST_DIR" 2>/dev/null)" ]; then \
        echo "[INFO] MLST scheme missing from build context, pulling from GitHub fork..."; \
        pseudoscope --pull-mlst-db; \
    else \
        echo "[OK] MLST scheme baked from build context"; \
        ls "$MLST_DIR" | head -3; \
    fi

RUN if [ ! -f /opt/pseudoscope/pseudoscope/modules/mge_module/mobileOG-db/beatrix-1-6_v1_all/mobileOG-db-beatrix-1.6.dmnd ]; then \
        echo "[INFO] mobileOG-db missing from build context, downloading..."; \
        pseudoscope --pull-mge-db; \
    fi

SHELL ["/bin/bash", "-c"]

RUN groupadd -g 1000 pseudoscope 2>/dev/null || true && \
    useradd -u 1000 -g 1000 -m -s /bin/bash pseudoscope 2>/dev/null || true && \
    chown -R 1000:1000 /opt/pseudoscope

ENV HOME=/tmp
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PSEUDOSCOPE_TMPDIR=/tmp

RUN printf '%s\n' \
    '#!/bin/bash' \
    'set -e' \
    'export CONDA_PREFIX=/opt/conda/envs/pseudoscope' \
    'export CONDA_DEFAULT_ENV=pseudoscope' \
    'export PATH=${CONDA_PREFIX}/bin:${PATH}' \
    'export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}' \
    'if [ $# -gt 0 ] && [ "${1#-}" != "$1" ]; then' \
    '    set -- pseudoscope "$@"' \
    'fi' \
    'exec "$@"' \
    > /usr/local/bin/entrypoint.sh \
    && chmod +x /usr/local/bin/entrypoint.sh

USER 1000:1000

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["pseudoscope", "-h"]
