FROM mambaorg/micromamba:0.14.0

RUN apt-get update && \
    apt-get install -y procps git bc && \
    rm -rf /var/lib/apt/lists/* && \
    CONDA_DIR="/opt/conda" && \
    git clone https://github.com/jimmyliu1326/InfluenzaNanopore.git && \
    chmod +x /InfluenzaNanopore/influenza_consensus.sh && \
    ln -s /InfluenzaNanopore/influenza_consensus.sh /usr/local/bin/influenza_consensus.sh

RUN micromamba install -n base -y -c bioconda -c conda-forge -f /InfluenzaNanopore/environment.yml && \
    micromamba clean --all --yes && \
    rm -rf $CONDA_DIR/conda-meta && \
    rm -rf $CONDA_DIR/include && \
    rm -rf $CONDA_DIR/lib/python3.*/site-packages/pip && \
    find $CONDA_DIR -name '__pycache__' -type d -exec rm -rf '{}' '+'

ADD db/ /db/