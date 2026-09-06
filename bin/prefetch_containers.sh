#!/usr/bin/env bash
#
# Pre-pull every Singularity container this pipeline uses, into the same
# shared cache Nextflow itself reads from on Gadi (singularity.cacheDir in
# conf/gadi.config). Run this ONCE from the Gadi LOGIN node (which has
# external network access) before ever launching the pipeline with
# -profile gadi — Gadi compute nodes have none, so a PBS-submitted task
# that needs to pull a container fresh will just fail.
#
# See docs/RUNNING_ON_GADI.md for the full picture.
#
# The output filenames below are NOT arbitrary — they reproduce Nextflow's
# own SingularityCache.simpleName() algorithm exactly (strip protocol,
# replace ':' and '/' with '-', add .img unless the source is itself a
# .sif), confirmed against Nextflow's source
# (SingularityCache.groovy) and cross-checked against a real example in a
# Nextflow GitHub discussion. If this drifts from a future Nextflow
# version's actual behaviour, Nextflow will just re-pull on its own
# (harmless on the login node, fatal on a compute node) — a mismatch is
# safe to notice because that re-pull attempt is exactly what this script
# exists to avoid.
#
# The four .sif entries already on Jayden's Gadi scratch (star/picard/gatk/
# samtools in nextflow.config) are plain local file paths, not URIs — they
# need no pulling and aren't listed here.
set -euo pipefail

PBS_PROJECT="${PBS_PROJECT:-jo11}"
CACHE_DIR="${SINGULARITY_CACHEDIR:-/scratch/${PBS_PROJECT}/${USER}/nxf_singularity_cache}"
mkdir -p "$CACHE_DIR"

echo "[INFO] Caching into: $CACHE_DIR"
echo "[INFO] Run this from the Gadi LOGIN node — compute nodes have no external network."
echo ""

# uri -> local cache filename, must match nextflow.config's containers{} block
declare -A IMAGES=(
    ["https://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_7"]="depot.galaxyproject.org-singularity-bwa-0.7.17--hed695b0_7.img"
    ["https://depot.galaxyproject.org/singularity/trim-galore:0.6.10--hdfd78af_0"]="depot.galaxyproject.org-singularity-trim-galore-0.6.10--hdfd78af_0.img"
    ["https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0"]="depot.galaxyproject.org-singularity-fastqc-0.12.1--hdfd78af_0.img"
    ["https://depot.galaxyproject.org/singularity/ensembl-vep:111.0--pl5321h2a3209d_0"]="depot.galaxyproject.org-singularity-ensembl-vep-111.0--pl5321h2a3209d_0.img"
    ["docker://griffithlab/pvactools:latest"]="griffithlab-pvactools-latest.img"
    ["https://depot.galaxyproject.org/singularity/star-fusion:1.13.0--hdfd78af_2"]="depot.galaxyproject.org-singularity-star-fusion-1.13.0--hdfd78af_2.img"
    ["https://depot.galaxyproject.org/singularity/arriba:2.5.1--h87b9561_0"]="depot.galaxyproject.org-singularity-arriba-2.5.1--h87b9561_0.img"
    ["https://depot.galaxyproject.org/singularity/telescope:1.0.3--py38h24c8ff8_1"]="depot.galaxyproject.org-singularity-telescope-1.0.3--py38h24c8ff8_1.img"
    ["docker://humanlongevity/hla:0.0.0"]="humanlongevity-hla-0.0.0.img"
)

for uri in "${!IMAGES[@]}"; do
    fname="${IMAGES[$uri]}"
    dest="${CACHE_DIR}/${fname}"
    if [ -f "$dest" ]; then
        echo "[SKIP] already cached: ${fname}"
        continue
    fi
    echo "[PULL] ${uri} -> ${fname}"
    singularity pull "$dest" "$uri"
done

echo ""
echo "[DONE] All containers cached. -profile gadi runs should now find them without a network call."
