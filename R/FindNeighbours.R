#' FindNeighbours
#' 
#' Masks Seurat FindNeighbors but with British English spelling
#'
#' @inheritParams Seurat::FindNeighbors
#' @inherit Seurat::FindNeighbors return
#' @export
#' 

FindNeighbours <- function(...){
  Seurat::FindNeighbors(...)
}