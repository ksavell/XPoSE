#' Save DimPlots for figures with specified colors
#'
#' @param seur_obj Seurat object
#' @param groupby Metadata column used for grouping
#' @param splitby Optional metadata column used for splitting
#' @param file_n Base file name
#' @param hex_list List of named hex color vectors
#' @param output_dir Directory where plots will be saved
#'
#' @return Saves PDF files
#' @export
#'
save_dimplot <- function(seur_obj,
                         groupby = NULL,
                         splitby = NULL,
                         file_n = NULL,
                         hex_list = NULL,
                         output_dir = ".") {

  # Ensure groupby and hex_list are provided
  if (is.null(groupby)) {
    stop("You must specify a groupby argument.")
  }

  if (is.null(hex_list) || !groupby %in% names(hex_list)) {
    stop("You must provide a valid hex list for the specified groupby.")
  }

  # Create output directory if needed
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # Get the hex color mapping for the specified groupby
  hex_colors <- hex_list[[groupby]]

  # Get split values, if requested
  split_values <- if (!is.null(splitby)) {
    unique(seur_obj@meta.data[[splitby]])
  } else {
    NULL
  }

  if (!is.null(split_values)) {

    for (split_val in split_values) {
      subset_obj <- seur_obj[
        ,
        seur_obj@meta.data[[splitby]] == split_val
      ]

      current_plot <- DimPlot(
        subset_obj,
        group.by = groupby,
        cols = hex_colors[
          names(hex_colors) %in%
            unique(subset_obj@meta.data[[groupby]])
        ],
        pt.size = 0.5,
        shuffle = TRUE
      ) +
        ggtitle(paste(groupby, "-", split_val)) +
        theme_void() +
        theme(legend.position = "none")

      output_file <- file.path(
        output_dir,
        paste0(
          file_n,
          "_",
          split_val,
          "_",
          groupby,
          ".pdf"
        )
      )

      pdf(
        file = output_file,
        width = 10,
        height = 10
      )

      print(current_plot)
      dev.off()
    }

  } else {

    current_plot <- DimPlot(
      seur_obj,
      group.by = groupby,
      cols = hex_colors[
        names(hex_colors) %in%
          unique(seur_obj@meta.data[[groupby]])
      ],
      pt.size = 0.5,
      shuffle = TRUE
    ) +
      ggtitle(groupby) +
      theme_void() +
      theme(legend.position = "none")

    output_file <- file.path(
      output_dir,
      paste0(
        file_n,
        "_",
        groupby,
        ".pdf"
      )
    )

    pdf(
      file = output_file,
      width = 10,
      height = 10
    )

    print(current_plot)

    dev.off()
  }
}
