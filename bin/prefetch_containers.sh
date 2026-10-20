#!/usr/bin/env bash
#
# Pre-pull every Singularity container this pipeline uses, into the same
# shared cache Nextflow itself reads from on Gadi (singularity.cacheDir in
# conf/gadi.config). Run this ONCE from the Gadi LOGIN node (which has
# external network access) before ever launching the pipeline with
# -profile gadi: Gadi compute nodes have none, so a PBS-submitted task
# that needs to pull a container fresh will just fail.
#
# See docs/RUNNING_ON_GADI.md for the full picture.
#
# The output filenames below are NOT arbitrary: they reproduce Nextflow's
# own SingularityCache.simpleName() algorithm exactly (strip protocol,
# replace ':' and '/' with '-', add .img unless the source is itself a
# .sif), confirmed against Nextflow's source
# (SingularityCache.groovy) and cross-checked against a real example in a
# Nextflow GitHub discussion. If this drifts from a future Nextflow
# version's actual behaviour, Nextflow will just re-pull on its own
# (harmless on the login node, fatal on a compute node): a mismatch is
# safe to notice because that re-pull attempt is exactly what this script
# exists to avoid.
#
# picard and gatk are plain local .sif paths on /scratch/jo11 and the
# HLA-HD image is built locally (containers/hlahd), so none of those are
# listed here.
#
# Usage:
#   bash bin/prefetch_containers.sh                      # pull what is missing
#   bash bin/prefetch_containers.sh --list-missing DIR   # print uncached images (used by run_somatic_gadi.sh)
set -euo pipefail

# must match nextflow.config's containers{} block
IMAGES=(
    "https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0"
    "https://depot.galaxyproject.org/singularity/trim-galore:0.6.10--hdfd78af_0"
    "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data"
    "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c2f02ceb8cbe5aa37930dcbc0de6f1d07dfb292178c7e780c44ab2f556bf829f/data"
    "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/164b5421d87e69f8e9f5add5a937e148e8165f11bfb18159da607f56ecdb1495/data"
    "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d837709891c2d98fc0956f6fd0dba18b0f67d96c4db74ccbae7db98fd00afe42/data"
    "https://depot.galaxyproject.org/singularity/ensembl-vep:111.0--pl5321h2a3209d_0"
    "docker://mgibio/bam_readcount_helper-cwl:1.1.1"
    "docker://griffithlab/vatools:5.1.1"
    "docker://griffithlab/pvactools:7.1.3"
)

# Nextflow's SingularityCache.simpleName(): drop the scheme, ':' and '/' -> '-', add .img
cache_name() { local n="${1#*://}"; n="${n//:/-}"; echo "${n//\//-}.img"; }

if [ "${1:-}" = "--list-missing" ]; then
    for uri in "${IMAGES[@]}"; do
        [ -f "${2:?cache dir}/$(cache_name "$uri")" ] || echo "$uri"
    done
    exit 0
fi

PBS_PROJECT="${PBS_PROJECT:-jo11}"
CACHE_DIR="${SINGULARITY_CACHEDIR:-/scratch/${PBS_PROJECT}/${USER}/nxf_singularity_cache}"
mkdir -p "$CACHE_DIR"

echo "[INFO] Caching into: $CACHE_DIR"
echo "[INFO] Run this from the Gadi LOGIN node: compute nodes have no external network."
echo ""

for uri in "${IMAGES[@]}"; do
    dest="${CACHE_DIR}/$(cache_name "$uri")"
    if [ -f "$dest" ]; then
        echo "[SKIP] already cached: $(basename "$dest")"
        continue
    fi
    echo "[PULL] ${uri}"
    singularity pull "$dest" "$uri"
done

echo ""
echo "[DONE] Public containers cached. HLA-HD still needs a local build: containers/hlahd/README.md"
