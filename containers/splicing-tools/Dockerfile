# Build from the repository root:
# docker build -t neoadjvax-splicing-tools:1.0 -f containers/splicing-tools/Dockerfile .
# SpliceAI code/model licensing must be reviewed before use or redistribution.
FROM mambaorg/micromamba:2.0.5
COPY --chown=$MAMBA_USER:$MAMBA_USER containers/splicing-tools/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN python -m pip install --no-cache-dir tensorflow==2.15.1 spliceai==1.3.1
# The regtools-aware Python 3 clusterer is pinned to an audited upstream commit.
ADD https://raw.githubusercontent.com/davidaknowles/leafcutter/2c9907ef66adf0bfb3092f0ceb6886ee5c046fbb/clustering/leafcutter_cluster_regtools.py /opt/leafcutter/leafcutter_cluster_regtools.py
# Nextflow overrides the image entrypoint. Make the environment available directly.
ENV PATH=/opt/conda/bin:$PATH
ENV TF_NUM_INTEROP_THREADS=1
RUN python -c "import spliceai; import tensorflow" && regtools --version && samtools --version && bcftools --version
