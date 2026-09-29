#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
JOBS_FILE="${1:-${SCRIPT_DIR}/cache_jobs.tsv}"
BUILD_SCRIPT="${SCRIPT_DIR}/00_build_nested_deg_cache.R"
VALIDATE_SCRIPT="${SCRIPT_DIR}/00_validate_cache_inputs.R"
LOG_DIR="${SCRIPT_DIR}/logs"
WORKER_FILE="${LOG_DIR}/.nested_cache_worker.sh"

CACHE_CPUS="${CACHE_CPUS:-8}"
CACHE_MEM="${CACHE_MEM:-48g}"
CACHE_TIME="${CACHE_TIME:-48:00:00}"
CACHE_LSCRATCH="${CACHE_LSCRATCH:-50}"
CACHE_MAX_CONCURRENT="${CACHE_MAX_CONCURRENT:-6}"
CACHE_CHUNK_SIZE="${CACHE_CHUNK_SIZE:-20}"

[[ -f "$JOBS_FILE" ]] || { echo "Jobs file not found: $JOBS_FILE" >&2; exit 2; }
[[ -f "$BUILD_SCRIPT" ]] || { echo "Build script not found: $BUILD_SCRIPT" >&2; exit 3; }
[[ -f "$VALIDATE_SCRIPT" ]] || { echo "Validation script not found: $VALIDATE_SCRIPT" >&2; exit 4; }

mkdir -p "$LOG_DIR"

if command -v module >/dev/null 2>&1; then
  module load R/4.5.2 || true
fi

Rscript "$VALIDATE_SCRIPT" "$JOBS_FILE"

N_JOBS=$(( $(wc -l < "$JOBS_FILE") - 1 ))
[[ "$N_JOBS" -eq 18 ]] || { echo "Expected 18 cache jobs; found $N_JOBS" >&2; exit 5; }

JOBS_ABS="$(cd "$(dirname "$JOBS_FILE")" && pwd -P)/$(basename "$JOBS_FILE")"
BUILD_ABS="$(cd "$(dirname "$BUILD_SCRIPT")" && pwd -P)/$(basename "$BUILD_SCRIPT")"

cat > "$WORKER_FILE" <<'WORKER'
#!/bin/bash
set -euo pipefail

JOBS_FILE="$1"
BUILD_SCRIPT="$2"
CHUNK_SIZE="$3"

if command -v module >/dev/null 2>&1; then
  module load R/4.5.2 || true
fi

line=$((SLURM_ARRAY_TASK_ID + 1))
region=$(awk -F '\t' -v n="$line" 'NR==n {print $1}' "$JOBS_FILE")
cluster=$(awk -F '\t' -v n="$line" 'NR==n {print $2}' "$JOBS_FILE")
pb_rds=$(awk -F '\t' -v n="$line" 'NR==n {print $3}' "$JOBS_FILE")
meta_rds=$(awk -F '\t' -v n="$line" 'NR==n {print $4}' "$JOBS_FILE")
cache_root=$(awk -F '\t' -v n="$line" 'NR==n {print $5}' "$JOBS_FILE")

[[ -n "$region" && -n "$cluster" && -n "$pb_rds" && -n "$meta_rds" && -n "$cache_root" ]] || {
  echo "Could not parse cache job row for array task $SLURM_ARRAY_TASK_ID" >&2
  exit 6
}

export TMPDIR="${LSCRATCH:-${TMPDIR:-/tmp}}"
export R_TEMP_DIR="$TMPDIR"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "region=$region cluster=$cluster"

Rscript "$BUILD_SCRIPT" \
  --region "$region" \
  --cluster "$cluster" \
  --pb_rds "$pb_rds" \
  --meta_rds "$meta_rds" \
  --cache_root "$cache_root" \
  --chunk_size "$CHUNK_SIZE" \
  --n_cores "${SLURM_CPUS_PER_TASK:-1}"
WORKER
chmod +x "$WORKER_FILE"

job_id=$(sbatch --parsable \
  --job-name=xpose_degcache \
  --array="1-${N_JOBS}%${CACHE_MAX_CONCURRENT}" \
  --cpus-per-task="$CACHE_CPUS" \
  --mem="$CACHE_MEM" \
  --time="$CACHE_TIME" \
  --gres="lscratch:${CACHE_LSCRATCH}" \
  --output="$LOG_DIR/cache-%A_%a.out" \
  --error="$LOG_DIR/cache-%A_%a.err" \
  "$WORKER_FILE" "$JOBS_ABS" "$BUILD_ABS" "$CACHE_CHUNK_SIZE")

echo "Submitted nested cache array: $job_id"
echo "Jobs: $JOBS_ABS"
echo "Array: 1-${N_JOBS}%${CACHE_MAX_CONCURRENT}"
