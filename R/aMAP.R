#' aMAP
#'
#' Creates an annotated UMAP plot with extended labels pointing to cluster centers.
#' This function adds lines and text labels extending from cluster centroids to white space, originating from the centre of the plot,
#' with customizable extension and text offset per cluster, and UMAP centre parameters.
#'
#' @param seurat_object A Seurat object containing UMAP reduction (Seurat::RunUMAP() must have already been run)
#' @param a_centre_1 X-coordinate of the plot center (default: calculated from data)
#' @param a_centre_2 Y-coordinate of the plot center (default: calculated from data)
#' @param default_ext Default extension factor for all clusters (default: 2)
#' @param per_exts Named numeric vector of extension factors per cluster (default: default_ext (NULL))
#' @param text_offset Default text offset from label positions (default: 0.5)
#' @param per_text_offsets Named numeric vector of text offsets per cluster (default: text_offset (NULL))
#' @param ... Additional parameters passed to Seurat::DimPlot() (All default Seurat::DimPlot() options are possible, including group.by, split.by, order, cols, pt.size)
#'
#' @return A ggplot2 object with annotated UMAP
#'
#' @examples
#' # EXAMPLE 1: Cluster ordering (for plotting order on UMAP, passed to Seurat::DimPlot())
#' clus_order <- rev(c(
#'   "NK",
#'   "CD4 T",
#'   "CD8 T",
#'   "Cycling T",
#'   "B",
#'   "Macrophage",
#'   "Monocyte",
#'   "DC",
#'   "Endothelial",
#'   "Lining Fibroblast",
#'   "Sublining Fibroblast",
#'   "Plasma",
#'   "Mural"
#' ))
#' # Basic usage:
#' # p1 <- aMAP(seurat_object = data, order = clus_order)
#' # p1
#'
#' # EXAMPLE 2: Custom extension factors per cluster (Not all clusters have to be defined)
#' # (Controls how far label lines extend from cluster centers)
#' per_exts <- c(
#'   "Macrophage" = 4,
#'   "Cycling T" = 1.75,
#'   "Monocyte" = 2.5,
#'   "Plasma" = 2.75,
#'   "B" = 1.5,
#'   "CD4 T" = 3.5,
#'   "CD8 T" = 4,
#'   "NK" = 2.5,
#'   "Endothelial" = 2.5,
#'   "Mural" = 1.5,
#'   "Sublining Fibroblast" = 4.25,
#'   "Lining Fibroblast" = 3.75,
#'   "DC" = 2
#' )
#' # With custom extension factors:
#' # p2 <- aMAP(seurat_object = data,
#' #            per_exts = per_exts)
#' # p2
#'
#' # EXAMPLE 3: Full customization with text offsets (Not all clusters have to be defined)
#' # (Controls distance from label line to text)
#' text_offsets <- c(
#'   "Macrophage" = 1.5,
#'   "Monocyte" = 2,
#'   "Plasma" = 1,
#'   "Endothelial" = 2.25
#' )
#' # Complete customization:
#' # p3 <- aMAP(seurat_object = data,
#' #            per_exts = per_exts,
#' #            per_text_offsets = text_offsets)
#' # p3
#'
#' @export

aMAP <- function(seurat_object, a_centre_1 = NULL, a_centre_2 = NULL,
                 default_ext = NULL, per_exts = NULL, text_offset = NULL,
                 per_text_offsets = NULL, ...){
  
  if(is.null(default_ext)){
    default_ext = 2
  }
  
  if(is.null(text_offset)){
    text_offset = 0.5
  }
  
  p <- DimPlot(..., object = seurat_object, label = F,
               label.box = F, raster = F, repel = F) + NoLegend() + ggtitle("") +
    NoAxes() + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
  
  # Extract the data from the plot
  umap_data <- p$data
  
  # Rename in case using a different key for the reduction
  names(umap_data)[1:2] <- c("umap_1", "umap_2")
  umap_data$umap_1 <- as.numeric(umap_data$umap_1)
  umap_data$umap_2 <- as.numeric(umap_data$umap_2)
  
  # Extract the grouping column name automatically
  group_name <- colnames(umap_data)[3]
  
  # Calculate centroids for each cluster
  formula <- as.formula(paste("cbind(umap_1, umap_2) ~", group_name))
  centroids <- aggregate(formula, data = umap_data, FUN = median)
  names(centroids)[2:3] <- c("center_umap_1", "center_umap_2")
  
  # Calculate the overall center of all points
  if(is.null(a_centre_1) | is.null(a_centre_2)){
    overall_center <- colMeans(umap_data[, c("umap_1", "umap_2")])
    print(paste0("a_centre_1 & a_centre_2 not specified, using calculated centre at ",
                 paste0(overall_center[1], ", ", overall_center[2])))
  } else {
    overall_center <- c(umap_1 = a_centre_1, umap_2 = a_centre_2)
  }
  
  # Calculate direction vectors from overall center to each centroid
  centroids$dir_umap_1 <- centroids$center_umap_1 - overall_center["umap_1"]
  centroids$dir_umap_2 <- centroids$center_umap_2 - overall_center["umap_2"]
  
  # Normalize the direction vectors
  dir_length <- sqrt(centroids$dir_umap_1^2 + centroids$dir_umap_2^2)
  
  # Handle cases where dir_length is 0 (centroid at overall_center)
  if (any(dir_length == 0)) {
    centroids$dir_umap_1[dir_length == 0] <- 1
    centroids$dir_umap_2[dir_length == 0] <- 0
    dir_length[dir_length == 0] <- 1
  }
  
  # Normalize direction vectors
  centroids$dir_umap_1 <- centroids$dir_umap_1 / dir_length
  centroids$dir_umap_2 <- centroids$dir_umap_2 / dir_length
  
  # Assign extension factors per cluster
  extension_factors <- per_exts
  centroids[[group_name]] <- as.character(centroids[[group_name]])
  
  if (is.null(extension_factors)) {
    print('Using default extension factor. Specify "per_exts" for per-cluster values.')
    centroids$extension_factor <- default_ext
  } else {
    centroids$extension_factor <- ifelse(centroids[, group_name] %in% names(extension_factors),
                                         extension_factors[centroids[, group_name]],
                                         default_ext)
  }
  
  # Calculate label positions
  centroids$label_umap_1 <- centroids$center_umap_1 + centroids$dir_umap_1 * centroids$extension_factor
  centroids$label_umap_2 <- centroids$center_umap_2 + centroids$dir_umap_2 * centroids$extension_factor
  
  # Assign text offsets per cluster
  if (is.null(per_text_offsets)) {
    print('Using default text offset. Specify "per_text_offsets" for per-cluster values.')
    centroids$text_offset <- text_offset
  } else {
    centroids$text_offset <- ifelse(centroids[, group_name] %in% names(per_text_offsets),
                                    per_text_offsets[centroids[, group_name]],
                                    text_offset)
  }
  
  # Calculate final text positions
  centroids$text_umap_1 <- centroids$label_umap_1 + centroids$dir_umap_1 * centroids$text_offset
  centroids$text_umap_2 <- centroids$label_umap_2 + centroids$dir_umap_2 * centroids$text_offset
  
  # Calculate plot limits to include all labels
  x_min <- min(c(p$data$umap_1, centroids$text_umap_1), na.rm = TRUE) - 1
  x_max <- max(c(p$data$umap_1, centroids$text_umap_1), na.rm = TRUE) + 1
  y_min <- min(c(p$data$umap_2, centroids$text_umap_2), na.rm = TRUE) - 1
  y_max <- max(c(p$data$umap_2, centroids$text_umap_2), na.rm = TRUE) + 1
  
  # Update plot with expanded limits and add annotations
  p <- p + xlim(x_min, x_max) + ylim(y_min, y_max) +
    geom_segment(data = centroids,
                 aes(x = center_umap_1, y = center_umap_2,
                     xend = label_umap_1, yend = label_umap_2),
                 color = "black", size = 1) +
    geom_text(data = centroids,
              aes(x = text_umap_1, y = text_umap_2, label = !!sym(group_name)),
              color = "black", size = 5.5, fontface = 1)
  
  return(p)
}

