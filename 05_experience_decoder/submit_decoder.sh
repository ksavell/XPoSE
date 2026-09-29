#!/bin/bash
set -euo pipefail

HPC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
MAIL_USER=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mail-user)
      MAIL_USER="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

export DEG_CACHE_ROOT="${DEG_CACHE_ROOT:-/data/$USER/hierarchical_paired_deg_cache}"
export CACHE_RUN_ID="${CACHE_RUN_ID:-paired_deg_cache_5pct}"
export TOP_FRAC="${TOP_FRAC:-0.05}"
export DECODER_OUT_ROOT="${DECODER_OUT_ROOT:-/data/$USER/xpose_experience_decoder}"
export DECODER_RUN_ID="${DECODER_RUN_ID:-experience_decoder}"
export COLLECTED_OUT_DIR="${COLLECTED_OUT_DIR:-$DECODER_OUT_ROOT/collected}"
export DECODER_PERM_CHUNK_SIZE="${DECODER_PERM_CHUNK_SIZE:-10}"

VM_CLUSTERS="ITL23,ITL5,ITL6,ITvm,CTL6,ETL5,NPL5,Pvalb,Sst,Vip"
DM_CLUSTERS="ITL23,ITL5,ITL6,CTL6,ETL5,NPL5,Pvalb,Sst"

mkdir -p "$DECODER_OUT_ROOT" "$HPC_DIR/logs"

module purge
module load R/4.5.2

Rscript -e 'files <- commandArgs(TRUE); for (f in files) {parse(file=f); cat("Parsed:", f, "\n")}' \
  "$HPC_DIR/01_run_decoder.R" \
  "$HPC_DIR/02_collect_decoder.R"

Rscript -e 'pkgs <- c("optparse","glmnet","dplyr","tidyr","readr","tibble"); miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; if(length(miss)) stop("Missing R packages: ", paste(miss, collapse=", "))'

frac_tag=$(Rscript -e 'x<-as.numeric(commandArgs(TRUE)[1]); p<-x*100; if(abs(p-round(p))<1e-10) cat(sprintf("top_%dpct",round(p))) else cat(paste0("top_",sub("0+$","",sub("\\.$","",sprintf("%.3f",p))),"pct"))' "$TOP_FRAC")

check_cache () {
  local region="$1"
  local cluster="$2"
  local file="$DEG_CACHE_ROOT/$CACHE_RUN_ID/$region/$cluster/$frac_tag/nested_deg_and_expression_cache.rds"
  if [[ ! -f "$file" ]]; then
    echo "ERROR: missing nested cache: $file" >&2
    return 1
  fi
}

IFS=',' read -r -a VM_ARRAY <<< "$VM_CLUSTERS"
IFS=',' read -r -a DM_ARRAY <<< "$DM_CLUSTERS"

for cluster in "${VM_ARRAY[@]}"; do
  check_cache vmPFC "$cluster"
done
for cluster in "${DM_ARRAY[@]}"; do
  check_cache dmPFC "$cluster"
done

echo "Cache preflight passed."

MANIFEST="$HPC_DIR/decoder_jobs.tsv"
printf "vmPFC\t%s\n" "$VM_CLUSTERS" > "$MANIFEST"
printf "dmPFC\t%s\n" "$DM_CLUSTERS" >> "$MANIFEST"

MAIL_ARGS=()
if [[ -n "$MAIL_USER" ]]; then
  MAIL_ARGS=(--mail-user="$MAIL_USER" --mail-type=END,FAIL)
fi

DECODER_JOB=$(sbatch --parsable "${MAIL_ARGS[@]}" \
  --array="1-2%2" \
  --output="$HPC_DIR/logs/decoder-%A_%a.out" \
  --error="$HPC_DIR/logs/decoder-%A_%a.err" \
  --export=ALL \
  "$HPC_DIR/decoder_worker.sbatch" "$HPC_DIR" "$MANIFEST")

echo "Decoder array: $DECODER_JOB"

COLLECT_JOB=$(sbatch --parsable "${MAIL_ARGS[@]}" \
  --dependency="afterok:${DECODER_JOB}" \
  --output="$HPC_DIR/logs/decoder-collect-%j.out" \
  --error="$HPC_DIR/logs/decoder-collect-%j.err" \
  --export=ALL \
  "$HPC_DIR/collect_worker.sbatch" "$HPC_DIR")

echo "Collector: $COLLECT_JOB"
echo "Collected results: $COLLECTED_OUT_DIR"
