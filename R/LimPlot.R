#' LimPlot
#'
#' Masks Seurat DimPlot but with better defaults
#'
#' @inheritParams Seurat::FindNeighbors
#' @inherit Seurat::FindNeighbors return
#' @export
#'

LimPlot <- function(...){
  Seurat::DimPlot(..., label = T,  label.box = T, repel = T) + NoLegend()
}
