#' avgHeatmap
#' Plots an averaged heatmap, different from the default Seurat DoHeatmap
#' Usage: this will show the heatmap, use ggplot2 to save this out
#'
#' @param seurat_obj A Seurat object to plot features for
#' @param features A list of features to plot
#' @param group.by A metadata column by which to split your data
#' @param assay The assay in the Seurat object to pull from
#' @return A ggplot2 scaled data heatmap
#' @export


avgHeatmap <- function(seurat_obj, features = NULL, group.by = NULL, assay = assay) {
  # Define default markers if none provided
  if (is.null(features)) {
    features <- c("CD3D", "CD4", "CD8A", "GZMB", "CD79A", "SDC1", "MS4A1", "CD68", "CLEC10A", "CLIC5", "PRG4", "PDPN", "THY1", "PECAM1", "CD34", "MKI67")  # Example immune markers; replace with your defaults
  }

  # Calculate average expression per group
  avg_exp <- AggregateExpression(object = seurat_obj, assays = assay, features = features, group.by = group.by)
  avg_exp <- avg_exp[[1]]
  avg_exp <- avg_exp[features, , drop = FALSE]

  # Scale for heatmap
  avg_exp_scaled <- t(scale(t(avg_exp)))

  #cluster_group <- unique(data$fin.clus)
  #names(cluster_group) <- colnames(avg_exp_scaled)
  annotation_col <- data.frame(Cluster = colnames(avg_exp_scaled))
  rownames(annotation_col) <- colnames(avg_exp_scaled)
  ann_vals <- unique(annotation_col$Cluster)
  ann_colors <- list(
    Cluster = setNames(
      c(brewer.pal(length(annotation_col$Cluster), "Set3") ),
      ann_vals
    )
  )

  # Plot heatmap
  p <<- pheatmap::pheatmap(
    avg_exp_scaled,
    border_color = NA,
    color = colorRampPalette(c("#FE11FF","#0B0801","#FFF957"))(8),
    show_rownames = TRUE,
    show_colnames = TRUE,
    cluster_rows = F,
    cluster_cols = F,
    fontsize = 25,
    fontsize_row = 25,
    scale = "none",
    angle_col = 45,
    cellwidth = 40,
    cellheight = 30,
    annotation_col = annotation_col,
    annotation_colors = ann_colors
  )
  return(p)
}
