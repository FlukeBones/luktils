#' LotPlot
#'
#' Creates a split up Dot plot for better visual clarity.
#' This function adds to Seurat's DotPlot by breaking genes into groups as specified by user.
#'
#' @param seurat_object A Seurat object
#' @param gene_groups Groups of genes to box together on the DotPlot, defaults to groups of 5 if not specified
#' @param group_size Number of genes per group when gene_groups is not specified. Defaults to 5.
#' @param ... Additional parameters passed to Seurat::DotPlot() (All default Seurat::DotPlot() options are possible, including group.by, split.by, order, cols, pt.size)
#'
#' @return A ggplot2 object with a split up DotPlot
#'
#' @examples
#' # EXAMPLE 1: Custom groupings format
#' # You can set the groups by giving numbers directly
#' # (Controls groups and the names given to those groups)
#' gene_groups <- list(
#'   "Proliferating" = 1:10,
#'   "Communication" = 11:14,
#'   "Memory" = 15:20)
#' # Usage:
#' # p <- LotPlot(seurat_object = data,
#' #            gene_groups = gene_groups)
#' # p
#'
#' # EXAMPLE 2: Custom groupings format
#' # You can also use the length() and unlist() to convert a list of concatenated characters automatically:
#' list_1 <- c("Gene_1", "Gene_2", "Gene_3")
#' list_2 <- c("Gene_4", "Gene_5", "Gene_6", "Gene_7")
#' list_3 <- c("Gene_8", "Gene_9", "Gene_10", "Gene_11")
#' gene_groups <- list("Proliferating" = 1:length(list_1),
#'                    "Communication" = length(list_1)+1:length(list_2),
#'                    "Memory" = length(list_1)+length(list_2)+1:length(list_3))
#' # Usage:
#' p <- LotPlot(data, features = unlist(c(list_1, list_2, list_3)), gene_groups = gene_groups)
#' p
#'
#' @export

LotPlot <- function(seurat_object, gene_groups = NULL, group_size = 5, ...){
  
  p <- DotPlot(seurat_object, ...) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  #p
  
  # Extract data from DotPlot
  plot_data <- p$data
  
  # PRESERVE ORIGINAL CELL TYPE ORDER:
  original_order <- unique(p$data$id)
  
  # Create gene groups based on ACTUAL gene count
  n_genes <- length(unique(plot_data$features.plot))
  genes_unique <- unique(plot_data$features.plot)
  
  # Set default groupings if not provided
  if(is.null(gene_groups)){
    n_groups <- ceiling(n_genes / group_size)
    
    gene_groups <- list()
    for(i in 1:n_groups){
      start_idx <- (i-1) * group_size + 1
      end_idx <- min(i * group_size, n_genes)
      gene_groups[[paste0("Group ", i)]] <- start_idx:end_idx
    }
  }
  
  # Set groupings:
  
  plot_data$gene_group <- NA
  
  # Assign groups using the list
  for (group_name in names(gene_groups)) {
    indices <- gene_groups[[group_name]]
    genes_in_group <- genes_unique[indices]
    plot_data$gene_group[plot_data$features.plot %in% genes_in_group] <- group_name
  }
  
  # Convert to ordered factor automatically using list order
  plot_data$gene_group <- factor(plot_data$gene_group, levels = names(gene_groups))
  
  # AUTOMATICALLY EXTRACT UNIQUE GROUPS IN ORDER THEY APPEAR:
  group_levels <- unique(plot_data$gene_group[!is.na(plot_data$gene_group)])
  
  # CONVERT TO ORDERED FACTOR:
  plot_data$gene_group <- factor(plot_data$gene_group, levels = group_levels)
  
  # Create faceted plot
  p <- ggplot(plot_data, aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled)) +
    geom_point() +
    facet_wrap(~gene_group, scales = "free_x", space = "free_x") +
    scale_size_continuous(range = c(1, 6), name = "% Expr") +
    scale_color_gradient2(low = "white", 
                          mid = "lightsteelblue1",
                          high = "midnightblue", 
                          midpoint = 0, name = "Avg Expr") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12, colour = "black", family = "sans"),
          axis.text.y = element_text(size = 12, colour = "black", family = "sans"),
          strip.text.x = element_text(size = 12, family = "sans"),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
          plot.background = element_rect(color = "black", fill = NA, linewidth = 1.5))
  p
  
  return(p)
}
