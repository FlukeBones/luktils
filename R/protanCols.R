#' protanCols
#' Protanomaly is red-green colourblindness characterised by red-weak vision. This gives a panel of colours that are better than Seurat default for colourblind viz.
#' Usage: cols <- protanCols assigns the colours to a new thing
#'
#' @return A list of colours
#' @export

# protanCols:
protanCols <- function(){
  list(
    "olivedrab",
    "aquamarine2",
    "darkkhaki",
    "tan4",
    "hotpink4",
    "rosybrown1",
    "orchid2",
    "purple2",
    "darkorange",
    "steelblue2",
    "indianred1",
    "mediumorchid2"
  )
}
