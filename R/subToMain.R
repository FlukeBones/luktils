#' subToMain
#' Move idents from subclustered Seurat objects to an object containing the same cells by matching the cell names
#' Usage: data_obj <- subToMain(sub = subclustered_obj, data_obj)
#'
#' @param sub A Seurat object subset from the main object with new cell idents to transfer
#' @param main A Seurat object containing the cell names present in the sub object
#' @return A Seurat object with subclustered idents in the Ident slot
#' @export

# subToMain:
subToMain <- function(sub, main){
  # Find common cell names between sub and main
  cc <- intersect(names(sub@active.ident), names(main@active.ident))

  # Combine levels from both active.ident factors
  cl <- union(levels(main@active.ident), levels(sub@active.ident[cc]))

  # Update the levels of main@active.ident
  levels(main@active.ident) <- cl

  # Update the identities in main for the common cells
  main@active.ident[cc] <- sub@active.ident[cc]

  return(main)
}
