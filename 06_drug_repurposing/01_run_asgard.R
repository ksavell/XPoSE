# ==============================================================================
# XPoSE-seq drug repurposing: vmPFC RT and NC Active vs Non-active
# ==============================================================================
# Runs Asgard for the two vmPFC activity contrasts used in the manuscript:
#   1. RT Active vs Non-active
#   2. NC Active vs Non-active
#
# Inputs are the final DESeq2 result tables from 03_main_de. All tested genes with
# finite log2 fold change are passed through the original mouse-symbol-to-human
# ortholog mapping workflow; no DESeq2 significance filter is applied before
# GetDrug. The custom therapeutic score retains the original ranked-signature
# reversal logic and cluster-proportion weighting.
#
# For RT, a compact selected-drug gene-level table is also exported for the four
# candidates used in the reversal-similarity heatmaps.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(homologene)
  library(Asgard)
  library(cmapR)
})

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

cfg <- list(
  de_root = "output/03_differential_expression/03_main_de",
  vm_meta_rds = "output/03_differential_expression/03_main_de/pseudobulk/main_meta_vmPFC.rds",
  output_root = "output/06_drug_repurposing/asgard",

  analyses = c(
    RT = "RT_active_RT_nonactive_vmPFC",
    NC = "NC_active_NC_nonactive_vmPFC"
  ),

  # Keep the same neuronal populations used by the final Asgard analysis.
  clusters = c(
    "ITL23", "ITL5", "ITL6", "ITvm", "CTL6", "CTL6b", "ETL5",
    "NPL5", "Pvalb", "Sst", "Sncg", "Vip", "Lamp5"
  ),

  # External Asgard/LINCS resources. These large files should remain outside
  # the repository; edit the paths here for the local installation.
  ref_location = "input/06_drug_repurposing/Asgard/DrugReference",
  ref_tissue = "central-nervous-system",
  lincs_tissue = "central nervous system",
  gse92742_gctx = "input/06_drug_repurposing/Asgard/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx",
  gse70138_gctx = "input/06_drug_repurposing/Asgard/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx",

  repurposing_unit = "drug",
  connectivity = "negative",
  drug_type = "FDA",
  fda_drugs_only = TRUE,
  min_treatments = 1,

  selected_drugs = c(
    "trazodone",
    "simvastatin",
    "noscapine",
    "dextromethorphan"
  ),

  # The similarity panels are RT-only, so the larger gene-level export is only
  # written for RT.
  gene_level_export_experience = "RT"
)

# ------------------------------------------------------------------------------
# Small utilities
# ------------------------------------------------------------------------------

combine_p_fisher <- function(p) {
  keep <- is.finite(p) & p > 0 & p <= 1
  if (sum(keep) < 2) return(NA_real_)
  stats::pchisq(
    -2 * sum(log(p[keep])),
    df = 2 * sum(keep),
    lower.tail = FALSE
  )
}

rank_percentile <- function(x) {
  n <- length(x)
  if (n <= 1) return(rep(0.5, n))
  (rank(-x, ties.method = "average") - 1) / (n - 1)
}

validate_file <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path)
  normalizePath(path)
}

# ------------------------------------------------------------------------------
# Read final DESeq2 results and convert to the Asgard human-gene query format
# ------------------------------------------------------------------------------

read_and_prepare_cluster_de <- function(analysis_name, cluster_name) {
  de_file <- file.path(
    cfg$de_root,
    analysis_name,
    "results",
    paste0(cluster_name, "_DESeq2_results.csv")
  )

  if (!file.exists(de_file)) return(NULL)

  de <- read_csv(de_file, show_col_types = FALSE) %>%
    transmute(
      mouse_gene = as.character(gene),
      score = as.numeric(log2FoldChange),
      P.Value = as.numeric(pvalue),
      adj.P.Val = as.numeric(padj)
    ) %>%
    filter(!is.na(score), is.finite(score), !is.na(mouse_gene), nzchar(mouse_gene))

  if (nrow(de) == 0) return(NULL)

  homologs <- homologene::homologene(
    genes = unique(de$mouse_gene),
    inTax = 10090,
    outTax = 9606
  )

  if (is.null(homologs) || nrow(homologs) == 0) return(NULL)

  homologs <- as.data.frame(homologs)[, 1:2, drop = FALSE]
  colnames(homologs) <- c("mouse_gene", "human_gene")

  mapped <- de %>%
    left_join(homologs, by = "mouse_gene") %>%
    filter(!is.na(human_gene), nzchar(human_gene)) %>%
    arrange(adj.P.Val) %>%
    distinct(human_gene, .keep_all = TRUE) %>%
    transmute(
      gene = human_gene,
      score,
      P.Value,
      adj.P.Val
    )

  if (nrow(mapped) == 0) return(NULL)

  mapped <- as.data.frame(mapped)
  rownames(mapped) <- mapped$gene
  mapped
}

# ------------------------------------------------------------------------------
# Cell-type proportions from the pseudobulk metadata
# ------------------------------------------------------------------------------

vm_meta <- readRDS(validate_file(cfg$vm_meta_rds, "vmPFC pseudobulk metadata"))

required_meta <- c("cluster_name", "n_nuclei")
missing_meta <- setdiff(required_meta, colnames(vm_meta))
if (length(missing_meta) > 0) {
  stop("vmPFC metadata missing columns: ", paste(missing_meta, collapse = ", "))
}

cluster_nuclei <- vm_meta %>%
  filter(cluster_name %in% cfg$clusters) %>%
  group_by(cluster_name) %>%
  summarise(n_nuclei = sum(as.numeric(n_nuclei), na.rm = TRUE), .groups = "drop")

if (nrow(cluster_nuclei) == 0 || sum(cluster_nuclei$n_nuclei) <= 0) {
  stop("No vmPFC nuclei were available for the configured Asgard clusters.")
}

cluster_prop <- setNames(
  100 * cluster_nuclei$n_nuclei / sum(cluster_nuclei$n_nuclei),
  cluster_nuclei$cluster_name
)

# ------------------------------------------------------------------------------
# Load shared Asgard/LINCS resources once
# ------------------------------------------------------------------------------

dir.create(cfg$output_root, recursive = TRUE, showWarnings = FALSE)

ref_dir <- validate_file(cfg$ref_location, "Asgard reference directory")
gctx_92742_path <- validate_file(cfg$gse92742_gctx, "GSE92742 GCTX")
gctx_70138_path <- validate_file(cfg$gse70138_gctx, "GSE70138 GCTX")

gene_info_file <- validate_file(
  file.path(ref_dir, paste0(cfg$ref_tissue, "_gene_info.txt")),
  "Asgard gene-info file"
)
drug_info_file <- validate_file(
  file.path(ref_dir, paste0(cfg$ref_tissue, "_drug_info.txt")),
  "Asgard drug-info file"
)
drug_response_file <- validate_file(
  file.path(ref_dir, paste0(cfg$ref_tissue, "_rankMatrix.txt")),
  "Asgard rank-matrix file"
)

message("Loading Asgard reference...")
my_gene_info <- read.table(gene_info_file, sep = "\t", header = TRUE, quote = "")
my_drug_info <- read.table(drug_info_file, sep = "\t", header = TRUE, quote = "")

drug_ref_profiles <- GetDrugRef(
  drug.response.path = drug_response_file,
  probe.to.genes = my_gene_info,
  drug.info = my_drug_info
)

asgard_reference_genes <- unique(rownames(drug_ref_profiles$drug.rank.matrix))

message("Loading LINCS GCTX matrices...")
gse92742 <- cmapR::parse_gctx(gctx_92742_path)
gse70138 <- cmapR::parse_gctx(gctx_70138_path)
gse92742_matrix <- gse92742@mat
gse70138_matrix <- gse70138@mat
rm(gse92742, gse70138)
invisible(gc())

data(cell_data, package = "Asgard")
data(col_meta_GSE92742, package = "Asgard")
data(col_meta_GSE70138, package = "Asgard")
data(gene_meta, package = "Asgard")
data(FDA.drug, package = "Asgard")

# ------------------------------------------------------------------------------
# Score one Asgard analysis and optionally export selected-drug gene-level data
# ------------------------------------------------------------------------------

score_asgard_run <- function(cluster_degs, cluster_drugs, export_gene_level = FALSE) {
  # One cluster-level row per drug, matching the original workflow.
  drug_rows <- list()

  for (cl in names(cluster_drugs)) {
    x <- cluster_drugs[[cl]]
    if (is.null(x) || nrow(x) == 0 || !cl %in% names(cluster_prop)) next

    x <- x[!duplicated(x$Drug.name), , drop = FALSE]
    if (isTRUE(cfg$fda_drugs_only)) {
      x <- x[x$Drug.name %in% FDA.drug, , drop = FALSE]
    }
    if (nrow(x) == 0) next

    drug_rows[[cl]] <- tibble(
      drug = as.character(x$Drug.name),
      cluster = cl,
      cluster_prop = unname(cluster_prop[cl]),
      p_value = as.numeric(x$P.value),
      fdr = as.numeric(x$FDR)
    )
  }

  if (length(drug_rows) == 0) stop("No valid Asgard drug rows were available.")
  drug_list <- bind_rows(drug_rows)

  candidate_drugs <- unique(drug_list$drug)

  # CNS LINCS signatures for candidate drugs.
  cell_lines <- subset(cell_data, primary_site == cfg$lincs_tissue)$cell_id
  meta_92742 <- subset(
    col_meta_GSE92742,
    cell_id %in% cell_lines & pert_iname %in% candidate_drugs
  )
  meta_70138 <- subset(
    col_meta_GSE70138,
    cell_id %in% cell_lines & pert_iname %in% candidate_drugs
  )

  common_cols <- intersect(colnames(meta_92742), colnames(meta_70138))
  drug_metadata <- rbind(
    meta_92742[, common_cols, drop = FALSE],
    meta_70138[, common_cols, drop = FALSE]
  )

  if (nrow(drug_metadata) == 0) stop("No CNS LINCS signatures matched the Asgard candidates.")

  common_gene_ids <- intersect(
    rownames(gse92742_matrix),
    rownames(gse70138_matrix)
  )

  sig_92742 <- intersect(drug_metadata$sig_id, colnames(gse92742_matrix))
  sig_70138 <- intersect(drug_metadata$sig_id, colnames(gse70138_matrix))

  resp_92742 <- gse92742_matrix[common_gene_ids, sig_92742, drop = FALSE]
  resp_70138 <- gse70138_matrix[common_gene_ids, sig_70138, drop = FALSE]
  response_matrix <- cbind(resp_92742, resp_70138)

  if (ncol(response_matrix) == 0) stop("No matching LINCS response signatures were found.")

  response_df <- as.data.frame(response_matrix, check.names = FALSE) %>%
    rownames_to_column("gene_id") %>%
    mutate(gene_id = as.character(gene_id))

  gene_map <- as.data.frame(gene_meta) %>%
    transmute(
      gene_id = as.character(pr_gene_id),
      human_gene = as.character(pr_gene_symbol)
    )

  response_df <- response_df %>%
    inner_join(gene_map, by = "gene_id") %>%
    filter(!is.na(human_gene), nzchar(human_gene))

  # Preserve the original one-row-per-symbol behavior after GCTX mapping.
  response_df <- as.data.frame(response_df)
  rownames(response_df) <- make.unique(response_df$human_gene)
  response_matrix_symbol <- as.matrix(
    response_df[, setdiff(colnames(response_df), c("gene_id", "human_gene")), drop = FALSE]
  )
  storage.mode(response_matrix_symbol) <- "numeric"

  # Treatment-count filter.
  treatment_counts <- tibble(
    drug = candidate_drugs,
    n_usable_sigs = vapply(
      candidate_drugs,
      function(drug) {
        sigs <- subset(drug_metadata, pert_iname == drug)$sig_id
        length(intersect(sigs, colnames(response_matrix_symbol)))
      },
      integer(1)
    )
  ) %>%
    mutate(pass = n_usable_sigs >= cfg$min_treatments)

  candidate_drugs <- treatment_counts$drug[treatment_counts$pass]
  drug_list <- drug_list %>% filter(drug %in% candidate_drugs)

  # Restrict the query to the same rank-matrix universe used by GetDrug.
  cluster_degs <- lapply(cluster_degs, function(x) {
    if (is.null(x) || nrow(x) == 0) return(x)
    keep <- intersect(rownames(x), asgard_reference_genes)
    x[keep, , drop = FALSE]
  })

  score_rows <- list()
  gene_level_rows <- list()

  for (drug in candidate_drugs) {
    treatments <- subset(drug_metadata, pert_iname == drug)$sig_id
    treatments <- intersect(treatments, colnames(response_matrix_symbol))
    if (length(treatments) < cfg$min_treatments) next

    mean_response <- rowMeans(
      response_matrix_symbol[, treatments, drop = FALSE],
      na.rm = TRUE
    )
    names(mean_response) <- rownames(response_matrix_symbol)

    drug_stats <- drug_list %>% filter(drug == !!drug)
    therapeutic_score <- 0

    for (cl in names(cluster_degs)) {
      cs <- drug_stats %>% filter(cluster == cl)
      query <- cluster_degs[[cl]]

      if (nrow(cs) == 0 || is.null(query) || nrow(query) == 0) next

      represented <- intersect(rownames(query), names(mean_response))
      if (length(represented) == 0) next

      neuronal_score <- as.numeric(query[represented, "score"])
      drug_response <- as.numeric(mean_response[represented])
      names(neuronal_score) <- represented
      names(drug_response) <- represented

      reversal_product <- -neuronal_score * drug_response
      reversed <- is.finite(reversal_product) & reversal_product > 0
      reversal_fraction <- sum(reversed) / length(represented)

      cluster_fdr <- as.numeric(cs$fdr[1])
      contribution <-
        (as.numeric(cs$cluster_prop[1]) / 100) *
        (-log10(pmax(cluster_fdr, .Machine$double.xmin))) *
        reversal_fraction

      therapeutic_score <- therapeutic_score + contribution

      if (isTRUE(export_gene_level) && tolower(drug) %in% cfg$selected_drugs) {
        neuronal_pct <- rank_percentile(neuronal_score)
        drug_pct <- rank_percentile(drug_response)
        displacement <- abs(neuronal_pct - drug_pct)

        gene_level_rows[[length(gene_level_rows) + 1]] <- tibble(
          drug = tolower(drug),
          cluster = cl,
          human_gene = represented,
          neuronal_log2FC = neuronal_score,
          mean_drug_response = drug_response,
          neuronal_rank = rank(-neuronal_score, ties.method = "average"),
          drug_rank = rank(-drug_response, ties.method = "average"),
          reversed = reversed,
          inverse_rank_strength = if_else(reversed, displacement, 0),
          neuronal_rank_percentile = neuronal_pct,
          drug_rank_percentile = drug_pct,
          rank_displacement = displacement,
          reversal_product = reversal_product,
          cluster_prop = as.numeric(cs$cluster_prop[1])
        )
      }
    }

    score_rows[[drug]] <- tibble(
      drug = drug,
      therapeutic_score = therapeutic_score
    )
  }

  score_tbl <- bind_rows(score_rows)

  combined_p <- tapply(drug_list$p_value, drug_list$drug, combine_p_fisher)
  combined_fdr <- p.adjust(combined_p, method = "BH")

  score_tbl <- score_tbl %>%
    mutate(
      p_value = unname(combined_p[drug]),
      fdr = unname(combined_fdr[drug])
    ) %>%
    arrange(desc(therapeutic_score))

  gene_level_tbl <- if (length(gene_level_rows) > 0) {
    bind_rows(gene_level_rows)
  } else {
    tibble()
  }

  list(
    scores = score_tbl,
    gene_level = gene_level_tbl,
    treatment_counts = treatment_counts
  )
}

# ------------------------------------------------------------------------------
# Run RT and NC only
# ------------------------------------------------------------------------------

run_summary <- list()

for (experience in names(cfg$analyses)) {
  analysis_name <- unname(cfg$analyses[experience])
  message("\n", strrep("=", 70))
  message("Running Asgard: ", experience, " | ", analysis_name)
  message(strrep("=", 70))

  out_dir <- file.path(cfg$output_root, experience)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cluster_degs <- list()
  for (cl in cfg$clusters) {
    prepared <- read_and_prepare_cluster_de(analysis_name, cl)
    if (!is.null(prepared) && nrow(prepared) > 0) {
      cluster_degs[[cl]] <- prepared
    }
  }

  if (length(cluster_degs) == 0) {
    stop("No usable DESeq2 result tables were found for ", analysis_name)
  }

  message("Clusters entering GetDrug: ", paste(names(cluster_degs), collapse = ", "))

  cluster_drugs <- GetDrug(
    gene.data = cluster_degs,
    drug.ref.profiles = drug_ref_profiles,
    repurposing.unit = cfg$repurposing_unit,
    connectivity = cfg$connectivity,
    drug.type = cfg$drug_type
  )

  scored <- score_asgard_run(
    cluster_degs = cluster_degs,
    cluster_drugs = cluster_drugs,
    export_gene_level = identical(experience, cfg$gene_level_export_experience)
  )

  write_csv(scored$scores, file.path(out_dir, "drug_scores.csv"))
  write_csv(scored$treatment_counts, file.path(out_dir, "drug_treatment_counts.csv"))

  if (nrow(scored$gene_level) > 0) {
    write_csv(
      scored$gene_level,
      file.path(out_dir, "selected_drug_ranked_signature_gene_level.csv.gz")
    )
  }

  run_summary[[experience]] <- tibble(
    experience = experience,
    analysis_name = analysis_name,
    n_clusters = length(cluster_degs),
    n_scored_drugs = nrow(scored$scores),
    n_selected_gene_level_rows = nrow(scored$gene_level)
  )
}

write_csv(bind_rows(run_summary), file.path(cfg$output_root, "asgard_run_summary.csv"))

message("\nDone. Asgard outputs: ", cfg$output_root)
