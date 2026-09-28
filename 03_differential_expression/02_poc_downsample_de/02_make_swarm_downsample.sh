#!/usr/bin/env bash
# XPoSE-seq: generate Biowulf swarm for POC downsampling DE analysis
#
# Usage:
#   bash 03_differential_expression/02_downsample_de/02_make_swarm_downsample.sh <run_id>
#   bash 03_differential_expression/02_downsample_de/02_make_swarm_downsample.sh <run_id> --submit
#
# Run from any location within a cloned XPoSE repository. The script resolves
# repository paths relative to its own location.

set -euo pipefail

RUN_ID="${1:-}"
SUBMIT_FLAG="${2:-}"

if [[ -z "$RUN_ID" || "$RUN_ID" == "--submit" ]]; then
  echo "ERROR: run_id required." >&2
  echo "Usage: bash 02_make_swarm_downsample.sh <run_id> [--submit]" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Analysis paths ----------------------------------------------------------
DATA_ROOT="${REPO_ROOT}/output/03_differential_expression/02_downsample_de/hpc_input"
OUT_ROOT="${REPO_ROOT}/output/03_differential_expression/02_downsample_de"
R_SCRIPT="${SCRIPT_DIR}/03_run_downsample_de_job.R"
DE_SCRIPT="${SCRIPT_DIR}/functions/single_factor_DESeq.R"
CLUSTERS_FILE="${SCRIPT_DIR}/clusters_kept.txt"

# Analysis settings -------------------------------------------------------
PERCENTAGES=(0 1 5 25 50 75 100)
ITERATIONS=100
ALPHA=0.05
MIN_CELL=1

# Biowulf resources -------------------------------------------------------
R_MODULE="R/4.5.2"
GB=16
THREADS=4
WALLTIME="24:00:00"
LSCRATCH="lscratch:10"

mapfile -t CLUSTERS < <(grep -v '^[[:space:]]*$' "$CLUSTERS_FILE")

SWARM_DIR="${OUT_ROOT}/swarm"
LOG_DIR="${OUT_ROOT}/swarm_logs"
mkdir -p "$SWARM_DIR" "$LOG_DIR"

SWARM="${SWARM_DIR}/downsample_de_${RUN_ID}.swarm"
: > "$SWARM"

for cl in "${CLUSTERS[@]}"; do
  for pct in "${PERCENTAGES[@]}"; do
    echo "Rscript ${R_SCRIPT} --cluster ${cl} --percentage ${pct} --data_root ${DATA_ROOT} --out_root ${OUT_ROOT} --de_script ${DE_SCRIPT} --run_id ${RUN_ID} --iterations ${ITERATIONS} --alpha ${ALPHA} --min_cell ${MIN_CELL}" >> "$SWARM"
  done
done

N=$(wc -l < "$SWARM")
echo "Wrote ${N} jobs to ${SWARM}"
echo "${#CLUSTERS[@]} cell types x ${#PERCENTAGES[@]} active fractions; ${ITERATIONS} iterations/job"

CMD=(
  swarm
  -f "$SWARM"
  -g "$GB"
  -t "$THREADS"
  --time "$WALLTIME"
  --gres "$LSCRATCH"
  --module "$R_MODULE"
  --job-name "ds_${RUN_ID}"
  --logdir "$LOG_DIR"
)

if [[ "$SUBMIT_FLAG" == "--submit" ]]; then
  "${CMD[@]}"
else
  echo
  echo "Submit with:"
  printf ' %q' "${CMD[@]}"
  echo
fi
