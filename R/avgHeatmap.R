#' @name avgHeatmap
#' @title avgHeatmap
#'
#' @description Plots an averaged heatmap, different from the default Seurat DoHeatmap
#'
#' @usage avgHeatmap(seurat_obj,features = NULL,group.by = NULL,assay = NULL,cols = NULL)
#'
#' @param seurat_obj A Seurat object to plot features for
#' @param features A list of features to plot
#' @param group.by A metadata column by which to split your data
#' @param assay The assay in the Seurat object to pull from
#' @param cols A vector of colors for the clusters (uses Seurat default if NULL)
#' @return A ggplot2 scaled data heatmap
#' @export

utils::globalVariables(c("colors", "p"))

avgHeatmap <- function(seurat_obj, features = NULL, group.by = NULL, assay = NULL, cols = NULL) {
  temp_seu <- seurat_obj
  temp_seu@meta.data$idents <- Seurat::Idents(temp_seu)

  # Define defaults
  if (is.null(features)) {
    features <- c("CD3D", "CD4", "CD8A", "GZMB", "CD79A", "SDC1", "MS4A1", "CD68", "CLEC10A", "CLIC5", "PRG4", "PDPN", "THY1", "PECAM1", "CD34", "MKI67")
  }
  if (is.null(group.by)) {
    group.by <- "idents"
  }

  if (is.null(assay)) {
    assay <- "RNA"
  }

  # Calculate average expression per group
  avg_exp <- Seurat::AggregateExpression(object = temp_seu, assays = assay, features = features, group.by = group.by)
  avg_exp <- as.matrix(avg_exp[[1]])
  avg_exp <- avg_exp[features, , drop = FALSE]

  # Scale for heatmap
  avg_exp_scaled <- t(scale(t(avg_exp)))

  #cluster_group <- unique(data$fin.clus)
  #names(cluster_group) <- colnames(avg_exp_scaled)
  annotation_col <- data.frame(Cluster = colnames(avg_exp_scaled))
  rownames(annotation_col) <- colnames(avg_exp_scaled)
  ann_vals <- unique(annotation_col$Cluster)

  # Generate or use provided colors
  if (is.null(cols)) {
    cols <- scales::hue_pal()(length(ann_vals))
  } else if (length(cols) < length(ann_vals)) {
    warning("Provided cols has fewer colors than clusters. Generating additional colors.")
    additional <- scales::hue_pal()(length(ann_vals) - length(cols))
    cols <- c(cols, additional)
  } else if (length(cols) > length(ann_vals)) {
    cols <- cols[1:length(ann_vals)]  # Trim excess
  }

  # Ensure names match exactly
  ann_colors <- list(
    Cluster = stats::setNames(cols, ann_vals)
  )

  # Validate colors
  if (!all(sapply(cols, function(x) is.character(x) && grepl("^#", x) || x %in% colors()))) {
    stop("Invalid colors in cols. Must be hex codes or named colors.")
  }

  # Plot heatmap
  p <- pheatmap::pheatmap(
    avg_exp_scaled,
    border_color = NA,
    color = grDevices::colorRampPalette(c("#FE11FF","#0B0801","#FFF957"))(8),
    show_rownames = TRUE,
    show_colnames = TRUE,
    cluster_rows = F,
    cluster_cols = F,
    fontsize = 12,
    fontsize_row = 12,
    scale = "none",
    angle_col = 45,
    cellwidth = 20,
    cellheight = 15,
    annotation_col = annotation_col,
    annotation_colors = ann_colors
  )
  return(p)
}
