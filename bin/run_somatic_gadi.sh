#!/bin/bash
#
# Launch the somatic neoantigen pipeline on Gadi with the Nextflow driver
# inside a detached tmux session. Every pipeline step
# is still its own PBS job (conf/gadi.config); tmux only keeps the
# lightweight driver alive after you log out. See docs/RUNNING_ON_GADI.md.
#
# Usage (from the repo root on Gadi, ideally inside a persistent session):
#   bash bin/run_somatic_gadi.sh my_params.yaml [extra nextflow args...]
#   bash bin/run_somatic_gadi.sh --check my_params.yaml   # preflight only
#
# Then: tmux attach -t neoadjvax_somatic   (Ctrl+B, D to detach again)
#
# Env overrides: PBS_PROJECT (jo11), TMUX_SESSION (neoadjvax_somatic),
# NXF_WORK (/scratch/$PBS_PROJECT/$USER/neoadjvax_work),
# NEXTFLOW_MODULE (nextflow), SINGULARITY_MODULE (singularity), NXF_PROFILE (gadi).
set -euo pipefail

CHECK_ONLY=false
if [ "${1:-}" = "--check" ]; then CHECK_ONLY=true; shift; fi
PARAMS_FILE="${1:?usage: bash bin/run_somatic_gadi.sh [--check] <params.yaml> [nextflow args...]}"
shift

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PARAMS_FILE="$(readlink -f "$PARAMS_FILE")"
PBS_PROJECT="${PBS_PROJECT:-jo11}"
SESSION="${TMUX_SESSION:-neoadjvax_somatic}"
PROFILE="${NXF_PROFILE:-gadi}"
export NXF_WORK="${NXF_WORK:-/scratch/${PBS_PROJECT}/${USER}/neoadjvax_work}"
# Keep the driver JVM small: it only submits and polls PBS jobs.
export NXF_OPTS="${NXF_OPTS:--Xms500m -Xmx3g}"
export NXF_ANSI_LOG=false
# Modules are not strict-syntax (v2) clean yet; Nextflow >= 25.10 defaults to v2.
export NXF_SYNTAX_PARSER=v1

fail=0
err()  { echo "[FAIL] $*"; fail=1; }
ok()   { echo "[ OK ] $*"; }
warn() { echo "[WARN] $*"; }

# --- modules -------------------------------------------------------------
if command -v module >/dev/null 2>&1; then
    module load "${NEXTFLOW_MODULE:-nextflow}" "${SINGULARITY_MODULE:-singularity}" \
        || err "module load nextflow/singularity failed (try NEXTFLOW_MODULE=nextflow/<version>)"
fi
command -v nextflow >/dev/null 2>&1 || err "nextflow not on PATH"
command -v tmux >/dev/null 2>&1 || err "tmux not on PATH"

# resourceLimits in conf/gadi.config needs Nextflow >= 24.04
if command -v nextflow >/dev/null 2>&1; then
    nxf_ver="$(nextflow -version 2>/dev/null | grep -oE 'version [0-9]+\.[0-9]+' | head -1 | cut -d' ' -f2)"
    major="${nxf_ver%%.*}"; minor="${nxf_ver#*.}"
    if [ -z "$nxf_ver" ]; then
        warn "could not read the Nextflow version"
    elif [ "$major" -gt 24 ] || { [ "$major" -eq 24 ] && [ "$minor" -ge 4 ]; }; then
        ok "Nextflow $nxf_ver"
    else
        err "Nextflow $nxf_ver is older than 24.04 (needed for resourceLimits); load a newer module"
    fi
fi

# --- params file ---------------------------------------------------------
[ -f "$PARAMS_FILE" ] || { echo "[FAIL] params file not found: $PARAMS_FILE"; exit 1; }
# `|| true`: a key absent from the file must give "", not trip set -e/pipefail
param() { { grep -E "^\s*$1\s*:" "$PARAMS_FILE" || true; } | head -1 | sed -E "s/^\s*$1\s*:\s*//; s/\s+#.*$//; s/^['\"]//; s/['\"]$//"; }

grep -q "CHANGE_ME" "$PARAMS_FILE" && err "params file still contains CHANGE_ME placeholders"

sheet="$(param samplesheet)"
if [ -z "$sheet" ]; then err "samplesheet not set"; elif [ -f "$sheet" ]; then ok "samplesheet: $sheet"; else err "samplesheet missing: $sheet"; fi
if [ -f "$sheet" ]; then
    n=0
    # every absolute path in the sheet (FASTQ/BAM/segments) must exist
    for f in $(tail -n +2 "$sheet" | tr ',' '\n' | tr -d '\r' | grep '^/' || true); do
        n=$((n + 1))
        [ -e "$f" ] || err "samplesheet references a missing file: $f"
    done
    ok "samplesheet: checked $n input paths"
    # rows without hla_alleles (last column) need the HLA-HD image
    if tail -n +2 "$sheet" | tr -d '\r' | grep -q ',$' && [ -z "$(param hlahd_container)" ]; then
        err "some rows have no hla_alleles but hlahd_container is not set (containers/hlahd/README.md)"
    fi
fi
hlahd="$(param hlahd_container)"
[ -z "$hlahd" ] || { [ -f "$hlahd" ] && ok "hlahd_container: $hlahd" || err "hlahd_container missing: $hlahd"; }

for key in vep_cache vep_plugins; do
    val="$(param "$key")"
    if [ -z "$val" ]; then err "$key not set"; elif [ -d "$val" ]; then ok "$key: $val"; else err "$key missing: $val"; fi
done
plugins="$(param vep_plugins)"
if [ -n "$plugins" ] && [ -d "$plugins" ]; then
    for pm in Wildtype.pm Frameshift.pm; do
        [ -f "$plugins/$pm" ] || err "vep_plugins has no $pm (pVACseq needs it)"
    done
fi
[ -n "$(param outdir)" ] || warn "outdir not set: results go to ./results"

# --- VCF resources: must exist with a .tbi beside them --------------------
vcf_ok() { if [ ! -f "$2" ]; then err "$1 missing: $2"; elif [ ! -f "$2.tbi" ]; then err "$1 has no index: $2.tbi"; else ok "$1: $2"; fi; }
val="$(param germline_resource)"; [ -n "$val" ] && vcf_ok germline_resource "$val" || err "germline_resource not set"
for key in panel_of_normals contamination_sites; do
    val="$(param "$key")"
    if [ -n "$val" ]; then vcf_ok "$key" "$val"; else warn "$key not set (recommended; see docs/DEPENDENCIES.md)"; fi
done

# --- reference files and local .sif images (nextflow.config defaults) -----
# A default only has to exist if the params file doesn't override it.
while read -r key default; do
    val="$(param "$key")"; val="${val:-$default}"
    case "$key" in
        known_*) vcf_ok "$key" "$val" ;;
        *) if [ -e "$val" ]; then ok "$key: $val"; else err "$key missing: $val (override it in the params file)"; fi ;;
    esac
done <<'DEFAULTS'
genome_fasta /scratch/jo11/GRCh38.fa/GRCh38.primary_assembly.genome.fa
known_dbsnp /scratch/jo11/NeoTrio_RNA/dbsnp_138.hg38.vcf.gz
known_mills /scratch/jo11/NeoTrio_RNA/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
known_1000g /scratch/jo11/NeoTrio_RNA/1000G_phase1.snps.high_confidence.hg38.vcf.gz
DEFAULTS
fasta="$(param genome_fasta)"; fasta="${fasta:-/scratch/jo11/GRCh38.fa/GRCh38.primary_assembly.genome.fa}"
[ -f "$fasta.bwt" ] || warn "no BWA index beside $fasta: the run will build one first (several CPU-hours)"
for f in picard.sif gatk_4.5.0.0.sif; do
    [ -f "/scratch/jo11/singularity_images/$f" ] || err "missing local image: /scratch/jo11/singularity_images/$f"
done

cache="${SINGULARITY_CACHEDIR:-/scratch/${PBS_PROJECT}/${USER}/nxf_singularity_cache}"
missing="$(bash "$REPO_DIR/bin/prefetch_containers.sh" --list-missing "$cache" || true)"
if [ -n "$missing" ]; then
    err "containers not cached in $cache (run: bash bin/prefetch_containers.sh):"
    echo "$missing" | sed 's/^/         /'
else
    ok "all public containers cached"
fi

case "$(hostname)" in
    gadi-login-*) warn "on a plain login node: NCI kills long-running login-node processes; prefer a persistent session (docs/RUNNING_ON_GADI.md)" ;;
esac

if [ "$fail" -ne 0 ]; then echo; echo "Preflight failed: fix the [FAIL] lines above."; exit 1; fi
echo; echo "Preflight passed."
$CHECK_ONLY && exit 0

# --- launch ----------------------------------------------------------------
if tmux has-session -t "$SESSION" 2>/dev/null; then
    echo "tmux session '$SESSION' already exists: tmux attach -t $SESSION (or set TMUX_SESSION)"
    exit 1
fi

mkdir -p "$NXF_WORK" "$REPO_DIR/logs"
log="$REPO_DIR/logs/somatic_$(date +%Y%m%d_%H%M%S).log"
extra=''; [ $# -eq 0 ] || printf -v extra ' %q' "$@"
# Exports and module loads go inside the command: tmux < 3.0 (Rocky 8)
# has no `new-session -e`, and the new shell does not inherit modules.
cmd="module load ${NEXTFLOW_MODULE:-nextflow} ${SINGULARITY_MODULE:-singularity}; export NXF_WORK=$(printf %q "$NXF_WORK") NXF_OPTS=$(printf %q "$NXF_OPTS") NXF_ANSI_LOG=false NXF_SYNTAX_PARSER=v1; cd $(printf %q "$REPO_DIR") && nextflow run main.nf -profile $(printf %q "$PROFILE") --pbs_project $(printf %q "$PBS_PROJECT") -params-file $(printf %q "$PARAMS_FILE") -work-dir $(printf %q "$NXF_WORK") -resume${extra} 2>&1 | tee $(printf %q "$log"); echo; echo \"Nextflow exited (log: $log). Press Enter to close.\"; read"

tmux new-session -d -s "$SESSION" bash -lc "$cmd"

echo "Started in tmux session '$SESSION'."
echo "  attach:  tmux attach -t $SESSION    (detach: Ctrl+B, then D)"
echo "  log:     $log"
echo "  jobs:    qstat -u $USER"
