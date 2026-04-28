#' @name freqCharts
#' @title freqCharts
#'
#' @description
#' Plots charts of any 2 metadata columns in a Seurat object on 3 different barcharts, split by raw numbers and frequency.
#'
#' @usage freqCharts(seurat_obj, meta.1 = NULL, meta.2 = NULL, cols = NULL)
#'
#' @param seurat_obj A Seurat object to plot abundance charts for
#' @param meta.1 First column in the metadata frame to use for X axis plotting, defaults to group column
#' @param meta.2 Second column in the metadata frame to use for Y axis plotting, defaults to active idents column
#' @param cols Colours to use to label meta.1
#' @return 2 ggplot2 chart objects showing raw numbers and scaled numbers
#' @import dplyr
#' @export

utils::globalVariables(c("p", "freq", "p1", "p2"))

freqCharts <- function(seurat_obj, meta.1 = NULL, meta.2 = NULL, cols = NULL) {
  temp_seu <- seurat_obj
  # Define default splits if none provided, will go by cluster and group column:
  if (is.null(meta.1)) {
    meta.1 <- "group"
  }
  if (is.null(meta.2)) {
    temp_seu@meta.data$idents <- Seurat::Idents(temp_seu)
    meta.2 <- "idents"
  }

  # Check if columns exist
  if (!meta.1 %in% colnames(temp_seu@meta.data)) {
    stop(paste("Column", meta.1, "not found in meta.data"))
  }
  if (!meta.2 %in% colnames(temp_seu@meta.data)) {
    stop(paste("Column", meta.2, "not found in meta.data"))
  }

  # Generate colors if not provided
  num_colors_needed <- length(unique(temp_seu@meta.data[[meta.2]]))
  if (is.null(cols)) {
    cols <- scales::hue_pal()(num_colors_needed)
  } else if (length(cols) < num_colors_needed) {
    warning("Provided cols has fewer colors than needed for meta.2 levels. Generating additional colors.")
    additional <- scales::hue_pal()(num_colors_needed - length(cols))
    cols <- c(cols, additional)
  }

  # Calculate proportions
  plot_df <- temp_seu@meta.data %>%
    dplyr::group_by(.data[[meta.1]], .data[[meta.2]]) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(.data[[meta.1]]) %>%
    dplyr::mutate(freq = n / sum(n)) %>%
    dplyr::ungroup()

  # Plot 1: Basic dodged bars
  #p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[meta.1]], y = freq, fill = .data[[meta.2]])) +
  #  ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
  #  ggplot2::scale_fill_manual(values = rev(cols)) +  # Added for consistent colors
  #  ggplot2::theme_bw() +
  #  ggplot2::labs(x = meta.2, y = "Proportion", fill = meta.1) +
  #  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[meta.1]], y = freq, fill = .data[[meta.2]])) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), color = "black") +  # Added black outlines
    ggplot2::scale_fill_manual(values = rev(cols)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1.0),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 12, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
      legend.title = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 14)
    ) +
    ggplot2::labs(x = meta.2, y = "Proportion", fill = meta.1) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  # Plot 2: Frequency with outlines
  p1 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[meta.1]], y = freq, fill = .data[[meta.2]])) +
    ggplot2::geom_col(color = "black") +
    ggplot2::scale_fill_manual(values = rev(cols)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1.0),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 12, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
      legend.title = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 14)
    ) +
    ggplot2::labs(x = "", y = "Frequency", fill = "Cluster") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 1, size = 16))

  # Plot 3: Absolute numbers
  options(scipen = 999)
  p2 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[meta.1]], y = n, fill = .data[[meta.2]])) +
    ggplot2::geom_col(color = "black") +
    ggplot2::scale_fill_manual(values = rev(cols)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1.0),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 12, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
      legend.title = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 14)
    ) +
    ggplot2::labs(x = "", y = "Count", fill = "Cluster") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 1, size = 16))

  p3 <- p/ p1 + p2
  return(p3)
}
