#' lUMAP
#' Different plotting for a UMAP, with points with an outline because it looks better
#' Usage: Makes UMAP which can be displayed by calling p, or saved with ggplot2
#'
#' @param seurat_obj A Seurat object to plot UMAP
#' @param cols A list of colours to pass to the UMAP, defaults to Seurat base colours if not specified
#' @return A ggplot2 UMAP
#' @export

lUMAP <- function(seurat_obj, cols = NULL){
  temp <- seurat_obj
  temp@meta.data$clusters <- Idents(temp)

  umap <- data.frame(rownames(temp@meta.data),Embeddings(temp, reduction = "umap")[,1], Embeddings(temp, reduction = "umap")[,2],
                     temp$clusters)
  colnames(umap) <- c("ID","UMAP1", "UMAP2", "cluster")
  # Calculate cluster centers
  centers <- umap %>%
    group_by(cluster) %>%
    summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

  # If cols is not specified, generate defaults
  if (is.null(cols)) {
    cols <- scales::hue_pal()(length(levels(Seurat::Idents(temp))))  # Cap at 24 colors
  }

  umap$color_col <- umap$cluster

  # Plot
  p <<- ggplot(umap, aes(x=UMAP1, y=UMAP2)) +
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
  return(p)
}
