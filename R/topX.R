#' @name topX
#' @title A simple solution to grab top X genes per cluster from a FindAllMarkers readout.
#' @usage topX(res, n_genes = 5, pct_cut = 0.1, p_cut = 0.05, cut_by = "pct.1", neg = F)
#'
#' @param res The FindAllMarkers readout that you want to select from.
#' @param n_genes Number of top genes per cluster to return. Default is 5.
#' @param pct_cut Minimum percent expression cutoff. Default is 0.1 (10%).
#' @param cut_by A column in the results dataframe by which to cut. Default is pct.1.
#' @param p_cut Adjusted p-value cutoff. Default is 0.05.
#' @param neg Binary which determines whether you want negative markers. Default is False.
#'
#' @return A tibble of top X genes for each cluster in your readout
#' @import dplyr
#' @export

utils::globalVariables(c("pct.1", "p_val_adj", "avg_log2FC", "cut_by"))

topX <- function(res, n_genes = 5, pct_cut = 0.1, p_cut = 0.05, cut_by = "pct.1", neg = F){
  if(neg == F){
    markers <- res %>%
      filter(.data[[cut_by]] >= pct_cut,
             p_val_adj < p_cut) %>%
      group_by(cluster) %>%
      slice_max(order_by = avg_log2FC, n = n_genes)
  }

  if(neg == T){
    markers <- res %>%
      filter(.data[[cut_by]] >= pct_cut,
             p_val_adj < p_cut) %>%
      group_by(cluster) %>%
      slice_max(order_by = -avg_log2FC, n = n_genes)
  }

  return(markers)
}
