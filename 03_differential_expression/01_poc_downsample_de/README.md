POC downsampling differential expression on Biowulf

This directory contains the Biowulf portion of the POC active-fraction downsampling analysis. The local preparation step splits the final POC Seurat object into one RDS file per cell population. Biowulf then runs one swarm job for each cell population x target active fraction, with 100 downsampling/DESeq2 iterations per job. After the swarm completes, the job outputs are combined on Biowulf and transferred back to the local machine for the final figure plot.

1. Prepare cluster-level input files locally

From the base XPoSE repository directory, run:

Rscript 03_differential_expression/01_poc_downsample_de/01_prep_split_clusters.R

This creates:

output/03_differential_expression/01_poc_downsample_de/hpc_input/

with one .rds file for each population listed in hpc/clusters_kept.txt.

2. Transfer the HPC inputs and scripts to Biowulf

Create a working directory on Biowulf, for example:

mkdir -p /data/$USER/xpose_downsample

Transfer the local hpc_input/ directory plus the following files from 03_differential_expression/01_poc_downsample_de/hpc/:

02_make_swarm_downsample.sh

03_run_downsample_de_job.R

03_single_factor_DESeq.R

04_collect_downsample.R

clusters_kept.txt

The Biowulf working directory should then look like:

xpose_downsample/
├── hpc_input/
│   ├── CTL6.rds
│   ├── ETL5.rds
│   ├── ITL23.rds
│   ├── ITL5.rds
│   ├── ITL6.rds
│   └── Sst.rds
├── 02_make_swarm_downsample.sh
├── 03_run_downsample_de_job.R
├── 03_single_factor_DESeq.R
├── 04_collect_downsample.R
└── clusters_kept.txt

3. Run the analysis on Biowulf

Log in to Biowulf and move to the working directory:

cd /data/$USER/xpose_downsample

The swarm uses R/4.5.2. The R environment must include the packages used by the analysis, including Seurat, optparse, readr, ggplot2, DESeq2, and the neurorestore/Libra package that provides to_pseudobulk().

First generate the swarm file and print the submission command without launching jobs:

bash 02_make_swarm_downsample.sh poc

This generates 42 swarm jobs: 6 cell populations x 7 target active fractions (0, 1, 5, 25, 50, 75, and 100%). Each job runs 100 downsampling/DESeq2 iterations. The current swarm settings request 16 GB memory, 4 threads, 24 hours, and 10 GB local scratch per job.

After checking the generated command, submit the swarm:

bash 02_make_swarm_downsample.sh poc --submit

Monitor the swarm with Biowulf job tools, for example:

sjobs

or:

squeue -u $USER

Job logs are written to:

swarm_logs/

and per-job analysis outputs are written under:

output/run_poc/

After all swarm jobs have completed, load R and combine the per-job results:

module load R/4.5.2
Rscript 04_collect_downsample.R --out_dir output/run_poc

The collected output includes de_summary_long.csv, de_tally_long.csv, and Prism-ready mean and standard-deviation tables.

4. Transfer the collected results back to the local machine

Transfer the completed Biowulf directory:

output/run_poc/

back to:

output/03_differential_expression/01_poc_downsample_de/run_poc/

in the local XPoSE repository.

Generate the final Figure 5 downsampling panel locally on macOS with:

Rscript 03_differential_expression/01_poc_downsample_de/05_plot_downsample.R \
  --out_dir output/03_differential_expression/01_poc_downsample_de/run_poc
