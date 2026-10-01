# ==============================================================================
# XPoSE-seq relapse GO analysis
# Generates UPREGULATED topGO tables used by 02_plot_go.R.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(topGO)
  library(AnnotationDbi)
})

if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  stop("Package 'org.Mm.eg.db' is required. Install with BiocManager::install('org.Mm.eg.db').")
}

# Configuration ----------------------------------------------------------
analyses <- list(
  RT_active_RT_nonactive_dmPFC = list(
    region = "dmPFC",
    results_dir = "output/03_differential_expression/03_main_de/RT_active_RT_nonactive_dmPFC/results"
  ),
  RT_active_RT_nonactive_vmPFC = list(
    region = "vmPFC",
    results_dir = "output/03_differential_expression/03_main_de/RT_active_RT_nonactive_vmPFC/results"
  )
)

output_root <- "output/04_relapse_transcriptional_characterization/03_go"
padj_threshold <- 0.05
go_background_baseMean_min <- 10
ontologies <- c("BP", "MF", "CC")
node_size <- 5
org_db_name <- "org.Mm.eg.db"

# Utilities --------------------------------------------------------------
parse_topgo_p <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^<[[:space:]]*", "", x)
  suppressWarnings(as.numeric(x))
}

org_db <- getExportedValue(org_db_name, org_db_name)
valid_symbols <- AnnotationDbi::keys(org_db, keytype = "SYMBOL")

map_symbols <- function(symbols) {
  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  symbols <- intersect(symbols, valid_symbols)

  if (length(symbols) == 0) {
    return(tibble(SYMBOL = character(), ENTREZID = character()))
  }

  AnnotationDbi::select(
    org_db,
    keys = symbols,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  ) %>%
    as_tibble() %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE)
}

run_cluster_go <- function(result_file, analysis_name, region, cluster) {
  res <- read_csv(result_file, show_col_types = FALSE)

  required <- c("gene", "baseMean", "log2FoldChange", "padj")
  missing <- setdiff(required, names(res))
  if (length(missing) > 0) {
    stop("Missing columns in ", result_file, ": ", paste(missing, collapse = ", "))
  }

  res <- res %>%
    transmute(
      gene = as.character(gene),
      baseMean = as.numeric(baseMean),
      log2FoldChange = as.numeric(log2FoldChange),
      padj = as.numeric(padj)
    )

  background_symbols <- res %>%
    filter(!is.na(baseMean), baseMean >= go_background_baseMean_min) %>%
    pull(gene) %>%
    unique()

  up_symbols <- res %>%
    filter(
      !is.na(baseMean), baseMean >= go_background_baseMean_min,
      !is.na(padj), padj < padj_threshold,
      log2FoldChange > 0
    ) %>%
    pull(gene) %>%
    unique()

  background_map <- map_symbols(background_symbols)
  up_map <- map_symbols(up_symbols)

  background_entrez <- unique(background_map$ENTREZID)
  up_entrez <- intersect(unique(up_map$ENTREZID), background_entrez)

  cluster_dir <- file.path(output_root, analysis_name, cluster, "02_GO")
  dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)

  if (length(background_entrez) == 0 || length(up_entrez) == 0) {
    return(tibble(
      analysis_name = analysis_name,
      region = region,
      cluster = cluster,
      ontology = NA_character_,
      n_background_symbols = length(background_symbols),
      n_background_entrez = length(background_entrez),
      n_up_degs = length(up_symbols),
      n_up_entrez = length(up_entrez),
      status = "skipped",
      reason = "No mapped GO background or upregulated DEG foreground"
    ))
  }

  gene_list <- factor(as.integer(background_entrez %in% up_entrez))
  names(gene_list) <- background_entrez

  entrez_to_symbol <- setNames(
    c(background_map$SYMBOL, up_map$SYMBOL),
    c(background_map$ENTREZID, up_map$ENTREZID)
  )
  entrez_to_symbol <- entrez_to_symbol[!duplicated(names(entrez_to_symbol))]

  summary_rows <- list()

  for (ontology in ontologies) {
    status <- "success"
    reason <- NA_character_

    tryCatch({
      go_data <- new(
        "topGOdata",
        ontology = ontology,
        allGenes = gene_list,
        geneSelectionFun = function(x) x == 1,
        annot = annFUN.org,
        mapping = org_db_name,
        ID = "entrez",
        nodeSize = node_size
      )

      n_nodes <- length(usedGO(go_data))
      if (numGenes(go_data) == 0 || numSigGenes(go_data) == 0 || n_nodes == 0) {
        status <- "skipped"
        reason <- "No usable ontology-specific genes or GO nodes"
      } else {
        elim <- runTest(go_data, algorithm = "elim", statistic = "fisher")
        classic <- runTest(go_data, algorithm = "classic", statistic = "fisher")

        members <- genesInTerm(go_data)
        sig_gene_map <- lapply(members, function(ids) {
          hits <- intersect(as.character(ids), up_entrez)
          syms <- unname(entrez_to_symbol[hits])
          syms <- sort(unique(syms[!is.na(syms) & nzchar(syms)]))
          if (length(syms) == 0) NA_character_ else paste(syms, collapse = ";")
        })

        sig_gene_tbl <- tibble(
          GO.ID = names(sig_gene_map),
          sig_genes = unname(unlist(sig_gene_map))
        )

        raw <- GenTable(
          go_data,
          elim = elim,
          classic = classic,
          topNodes = n_nodes
        ) %>%
          as_tibble() %>%
          rename(p_elim = elim, p_classic = classic) %>%
          mutate(
            p_elim = parse_topgo_p(p_elim),
            p_classic = parse_topgo_p(p_classic),
            Annotated = as.numeric(Annotated),
            Significant = as.numeric(Significant),
            Expected = as.numeric(Expected),
            EnrichmentScore = if_else(Expected > 0, Significant / Expected, NA_real_),
            BH_elim = p.adjust(p_elim, method = "BH"),
            BH_classic = p.adjust(p_classic, method = "BH"),
            cluster = cluster,
            region = region,
            direction = "up"
          ) %>%
          left_join(sig_gene_tbl, by = "GO.ID")

        write_csv(
          raw,
          file.path(cluster_dir, paste0(cluster, "_", ontology, "_up_combined_raw.csv"))
        )
      }
    }, error = function(e) {
      status <<- "failed"
      reason <<- conditionMessage(e)
    })

    summary_rows[[ontology]] <- tibble(
      analysis_name = analysis_name,
      region = region,
      cluster = cluster,
      ontology = ontology,
      n_background_symbols = length(background_symbols),
      n_background_entrez = length(background_entrez),
      n_up_degs = length(up_symbols),
      n_up_entrez = length(up_entrez),
      status = status,
      reason = reason
    )
  }

  bind_rows(summary_rows)
}

# Run both relapse analyses ---------------------------------------------
run_summary <- list()

for (analysis_name in names(analyses)) {
  info <- analyses[[analysis_name]]
  files <- list.files(
    info$results_dir,
    pattern = "_DESeq2_results\\.csv$",
    full.names = TRUE
  )

  if (length(files) == 0) {
    warning("No DE result files found for ", analysis_name)
    next
  }

  for (file in files) {
    cluster <- sub("_DESeq2_results\\.csv$", "", basename(file))
    message(analysis_name, " | ", cluster)

    run_summary[[paste(analysis_name, cluster, sep = "__")]] <- run_cluster_go(
      result_file = file,
      analysis_name = analysis_name,
      region = info$region,
      cluster = cluster
    )
  }
}

summary_tbl <- bind_rows(run_summary)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
write_csv(summary_tbl, file.path(output_root, "GO_run_summary.csv"))

message("GO analysis complete: ", output_root)
