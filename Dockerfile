FROM condaforge/mambaforge:4.9.2-5

RUN apt-get update && \
    apt-get install -y procps git bc && \
    rm -rf /var/lib/apt/lists/* && \
    CONDA_DIR="/opt/conda" && \
    git clone https://github.com/jimmyliu1326/InfluenzaNanopore.git && \
    chmod +x /InfluenzaNanopore/influenza_consensus.sh && \
    ln -s /InfluenzaNanopore/influenza_consensus.sh /usr/local/bin/influenza_consensus.sh

RUN mamba create -p $CONDA_DIR/envs/medaka -c conda-forge -c bioconda -y medaka --strict-channel-priority && \
    mamba env update -n base -f /InfluenzaNanopore/environment.yml && \
    mamba clean --all --yes && \
    rm -rf $CONDA_DIR/conda-meta && \
    rm -rf $CONDA_DIR/include && \
    rm -rf $CONDA_DIR/lib/python3.*/site-packages/pip && \
    find $CONDA_DIR -name '__pycache__' -type d -exec rm -rf '{}' '+' && \
    conda run -n medaka medaka tools download_models
    
RUN mkdir /data && \
    chmod a+rwX /data

ADD db/ /db/

WORKDIR /data

ENV XDG_CACHE_HOME=/data