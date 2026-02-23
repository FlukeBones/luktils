#' synPanel
#' A list of genes for
#' Usage: feats <- synPanel(panel) assigns a list of relevant synovial genes to your object
#'
#' @param panel The panel you want to call
#' @return A list of relevant genes for synovial tissue clustering
#' @export

synPanel <- function(panel = NULL){
  if(is.null(panel)){
    #default panel:
    list(
    "CD3D",
    "CD4",
    "CD8A",
    "CD68",
    "THY1",
    "PRG4",
    "CD79A",
    "MS4A1",
    "SDC1",
    "IL3RA",
    "GNLY",
    "PECAM1"
    )} else
   if(panel == "fibroblast" | panel == "Fibroblast" | panel == "fib" | panel == "Fib"){
     list(
     "PDPN", #
     "THY1", # Sublining
     "PRG4", # Lining
     "CLIC5", # Lining
     "MMP3", # Lining
     "IGFBP5", # Lining
     "MFAP5", # Sublining F2 CD34+
     "CD34", # Sublining F2 CD34+
     "PI16", # Sublining F2 CD34+
     "ACKR3", # Sublining F2 CD34+
     "DKK3", # Sublining F3/F4 DKK3+
     "CXCL12", # Sublining F4 CXCL12+
     "NOTCH3" # Sublining F7 NOTCH3+
     )} else
  if(panel == "CD8" | panel == "eight" | panel == "eights" | panel == "CD8s"){
    list(
      "CD3D",
      "CD3E",
      "CD4",
      "CD8A",
      "CD8B",
      "GZMA",
      "GZMB",
      "GZMH",
      "GZMK",
      "PRF1",
      "GNLY",
      "IFNG",
      "IL7R",
      "CCR7",
      "AREG",
      "SELL",
      "CD69",
      "PRDM1",
      "IKZF1",
      "CCL5",
      "TRDC",
      "MKI67",
      "STMN1"
    )} else
  if(panel == "CD4" | panel == "four" | panel == "fours" | panel == "CD4s"){
    list(
      "CD3D",
      "CD3E",
      "CD4",
      "CD8A",
      "GZMA",
      "GZMB",
      "GZMH",
      "GZMK",
      "PDCD1",
      "IL7R",
      "CD69",
      "LAG3",
      "PRDM1",
      "CXCL13",
      "FOXP3",
      "IL2RA",
      "TIGIT",
      "IKZF2",
      "CCR7",
      "SELL",
      "AREG",
      "MKI67",
      "STMN1",
      "IKZF1",
      "KLF12"
    )} else
  if(panel == "endo" | panel == "endothelial" | panel == "Endo" | panel == "Endothelial"){
    list(
      "PECAM1",
      "SPARC",
      "COL4A1",
      "LIFR",
      "ICAM1",
      "NOTCH4",
      "CXCL12",
      "LYVE1"
    )} else
  if(panel == "B" | panel == "bcell" | panel == "BCell" | panel == "Bcell" | panel == "plasma" | panel == "Plasma"| panel == "plasma cell" | panel == "Plasma cell"){
    list(
      "CD79A",
      "MS4A1",
      "CD52",
      "CD69",
      "CXCR4",
      "LTB",
      "CD52",
      "TCL1A",
      "ITGAX",
      "LAMP1",
      "SDC1",
      "JCHAIN",
      "IGHM",
      "IGKC",
      "IGHG3",
      "IGHG1"
    )} else
  if(panel == "macrophage" | panel == "myeloid" | panel == "DC" | panel == "Macro" | panel == "Myeloid" | panel == "Mye" | panel == "Mono"){
    list(
      "CD68",
      "MERTK",
      "SELENOP",
      "LYVE1",
      "MAF",
      "C1QA",
      "MRC1",
      "HBEGF",
      "SPP1",
      "C15orf48",
      "LGALS1",
      "LDHA",
      "MARCO",
      "VSIG4",
      "ISG15",
      "IL1B",
      "FCN1",
      "PLCG2",
      "CLEC9A",
      "IRF8",
      "CLEC10A",
      "LYZ",
      "FCER1A",
      "CD1C",
      "IFITM2",
      "CD48",
      "LY6E",
      "LAMP3",
      "IL7R",
      "CCR7",
      "IL3RA",
      "IRF7",
      "IRF8",
      "JCHAIN",
      "GZMB"
    )}
}

feats <- synPanel("macrophage")
head(feats)
