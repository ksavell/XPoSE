# 06_drug_repurposing

## Scripts

1. `01_run_asgard.R` reads the final DESeq2 result tables from `03_differential_expression/03_main_de`, maps the query genes to human orthologs, runs Asgard for RT and NC, computes the final cross-cell-type therapeutic scores/FDR values, and exports RT gene-level ranked-signature values for trazodone, simvastatin, noscapine, and dextromethorphan.
2. `02_plot_reversal_candidates.R` compares NC and RT final Asgard FDR values and produces the candidate scatter used in the manuscript.
3. `03_plot_reversal_similarity.R` uses the RT selected-drug gene-level export to calculate binary reversed-gene Jaccard and inverse-rank-strength weighted Jaccard similarity and produces the two selected-drug heatmaps.

## Required external Asgard resources

The LINCS GCTX matrices and the Asgard drug-reference files are large external resources and are not included in the repository. Edit the paths in `01_run_asgard.R` if these files are stored elsewhere.

Expected reference files:

- `central-nervous-system_gene_info.txt`
- `central-nervous-system_drug_info.txt`
- `central-nervous-system_rankMatrix.txt`
- `GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx`
- `GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx`

Run scripts in numerical order.
