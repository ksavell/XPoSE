# XPoSE-seq experience decoder

## Workflow

1. `00_make_cache_jobs.R` creates the 18 nested-cache jobs: 10 vmPFC cell types and 8 dmPFC cell types.
2. `00_validate_cache_inputs.R` verifies the pseudobulk/meta inputs and the complete 12-rat RT/NC Active/Non-active structure before expensive computation.
3. `00_submit_cache.sh` submits the 18 Biowulf cache jobs. Each job runs `00_build_nested_deg_cache.R` for one region x cell type.
4. `01_run_decoder.R` runs the final nested activity-difference decoder separately for dmPFC and vmPFC.
5. `02_collect_decoder.R` combines the two regional decoder summaries.
6. `03_plot_decoder_accuracy.R` plots percent correct for dmPFC and vmPFC.

## Cache inputs

The cache builder reads the raw pseudobulk counts and metadata produced by:

`03_differential_expression/03_main_de/01_create_pseudobulk_meta.R`

Expected default inputs are:

- `output/03_differential_expression/03_main_de/pseudobulk/main_pb_dmPFC.rds`
- `output/03_differential_expression/03_main_de/pseudobulk/main_meta_dmPFC.rds`
- `output/03_differential_expression/03_main_de/pseudobulk/main_pb_vmPFC.rds`
- `output/03_differential_expression/03_main_de/pseudobulk/main_meta_vmPFC.rds`

## Nested cache design

For each region x cell type, the cache builder enumerates all 600 sex-balanced five-rat training subsets that can arise from the same-sex paired holdout / exact label-assignment design. For each subset, Active versus Non-active DEGs are recomputed from raw pseudobulk counts with DESeq2 using the paired design `~ decoder_cond + decoder_rat`, with BH-adjusted P < 0.05 and no additional log2FC threshold.

The builder also creates the 30 possible same-sex held-out-pair expression caches. Activity-difference expression is normalized with a training-reference median-ratio procedure, filtered using training samples only, transformed as log2(normalized count + 1), differenced within rat as Active minus Non-active, adjusted for sex using coefficients estimated from the training rats only, and ordered by training-set variance.

The resulting cache file is:

`<DEG_CACHE_ROOT>/paired_deg_cache_5pct/<region>/<cluster>/top_5pct/nested_deg_and_expression_cache.rds`

Its decoder-facing components are `rat_meta`, `deg_by_subset`, and `expression_pair_cache`.

## Build caches on Biowulf

From the repository root:

```bash
export DEG_CACHE_ROOT=/data/$USER/hierarchical_paired_deg_cache
Rscript 05_experience_decoder/00_make_cache_jobs.R
Rscript 05_experience_decoder/00_validate_cache_inputs.R
bash 05_experience_decoder/00_submit_cache.sh
```

Default cache resources are 8 CPUs, 48 GB RAM, 48 hours, 50 GB local scratch, and a maximum of 6 simultaneous array jobs. Override with `CACHE_CPUS`, `CACHE_MEM`, `CACHE_TIME`, `CACHE_LSCRATCH`, `CACHE_MAX_CONCURRENT`, or `CACHE_CHUNK_SIZE`.

Cache construction is resumable at the DESeq2 chunk level. If a completed final cache already exists, the builder skips it unless `--rebuild` is supplied directly to `00_build_nested_deg_cache.R`.

## Final decoder

The final decoder uses only the `activity_difference` representation. Within each outer fold and exact RT/NC label assignment, it retrieves the Active-vs-Non-active DEG sets for the five rats assigned to the RT training class and the five rats assigned to the NC training class, takes their regional RT union NC union across configured cell types, and fits the nested hierarchical ridge decoder using only training-derived features and preprocessing.

The exact assignment space contains 400 sex-balanced RT/NC label assignments. The observed assignment contains 18 eligible same-sex RT/NC held-out pairs. Percent correct is `100 x pair concordance`, where a held-out pair is correct when the RT rat receives a higher RT probability than the paired NC rat.
