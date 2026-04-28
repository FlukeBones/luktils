#' @name lUMAP
#' @title lUMAP
#'
#' @description
#' Different plotting for a UMAP, with points with an outlin because it looks better for low #s of cells.
#'
#' @usage lUMAP(seurat_obj, cols = NULL)
#'
#' @param seurat_obj A Seurat object to plot UMAP
#' @param cols A list of colours to pass to the UMAP, defaults to Seurat base colours if not specified
#' @return A ggplot2 UMAP
#' @import dplyr
#' @import ggplot2
#' @export

utils::globalVariables(c("cluster", "UMAP1", "UMAP2", "p", "color_col"))

lUMAP <- function(seurat_obj, cols = NULL){
  temp <- seurat_obj
  temp@meta.data$clusters <- Seurat::Idents(temp)

  umap <- data.frame(rownames(temp@meta.data), Seurat::Embeddings(temp, reduction = "umap")[,1], Seurat::Embeddings(temp, reduction = "umap")[,2],
                     temp$clusters)
  colnames(umap) <- c("ID","UMAP1", "UMAP2", "cluster")
  # Calculate cluster centers
  centers <- umap %>%
    group_by(cluster) %>%
    dplyr::summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

  # If cols is not specified, generate defaults
  if (is.null(cols)) {
    cols <- scales::hue_pal()(length(levels(Seurat::Idents(temp))))
  }

  umap$color_col <- umap$cluster

  # Plot
  p <- ggplot(umap, aes(x=UMAP1, y=UMAP2)) +
    geom_point(aes(fill=factor(color_col)), color="black", size=2.5, shape=21, stroke=0.2) +
    scale_fill_manual(
      name = NULL,
      labels = levels(factor(umap$color_col)),
      values = cols
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    guides(fill = guide_legend(reverse = T, override.aes = list(shape = 22, size = 10)))+
    theme(legend.text = element_text(size = 14),legend.title = element_text(size = 14))
  rm(temp)
  return(p)
}
