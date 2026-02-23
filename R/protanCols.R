#' protanCols
#'
#' Protanomaly is red-green colourblindness characterised by red-weak vision. This gives a panel of colours that are better than Seurat default for colourblind viz.
#' Usage: cols <- protanCols assigns the colours to a new thing
#'
#' @param panel Panel to return, one of NULL (blank) "default" or "0", "1" or "2"
#' @return A list of colours
#' @export

# protanCols:
protanCols <- function(panel = NULL){
  if(is.null(panel) || panel == "default" || panel == "0"){
    c(
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
  )} else
    if(panel == "1"){
      c(
        "darkorange",
        "orchid1",
        "violetred1",
        "indianred1",
        "gold2",
        "firebrick4",
        "gold4"
        )} else
    if(panel == "2"){
      c(
        "#5C2175",
        "#D123CE",
        "#9447D3",
        "#6E9CBB",
        "dodgerblue3",
        "purple3",
        "cyan4"
        )} else
          stop('No such panel. Choose one of NULL / "default" / "0", "1" or "2".')
}
