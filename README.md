## luktils (luke-tils, like Luke Utils) is a package made to make my life easier.

<sub><sup>Yes that's how you say it. It seemed like a good idea at the time.</sub></sup>

There's no guarantee it'll make your life easier.

<sub><sup>Or mine, for that matter.</sub></sup>

scRNA-seq analysis, plotting, and dataviz.

<sub><sup>And it doesn't even do that particularly well.</sub></sup>


## Current Functions:
### / [subToMain](R/subToMain.R) \\
This is a quick function to move the active idents from an object subclustered new labels to a larger object by matching cell names.

### / [avgHeatmap](R/avgHeatmap.R) \\
This is a heatmap plotting that uses pheatmap and AggregateExpression to make averaged expression heatmaps per cluster, instead of the default DoHeatmap which has each cell plotted as it's own line. This uses the same heatmap colours as default Seurat plotting for clarity, colours specify the clustering bar along the top of the heatmap.

The following plot can be made by running:
> p <- avgHeatmap(data, features = synPanel("eight"), cols = protanCols(2))
> p
![luktils::avgHeatmap](images/avgHeatmap.png)

### / [freqCharts](R/freqCharts.R) \\
This makes 3 charts which show one metadata against another, both in raw numbers and in frequency (scaled from 0-1) to show, for example, abundance of clusters per experimental group.

### / [aMAP](R/aMAP.R) \\
aMAP (Annotated MAP) builds on Seurat's DimPlot to instead draw lines to annotate clusters, resulting in a cleaner UMAP without overlayed cluster names.

### / [lUMAP](R/lUMAP.R) \\
lUMAP (El-UMAP or Loo-MAP, as in "Luke UMAP") uses a previously generated UMAP graph but changes the plot size, point circling, larger legend, and plot boxed in.

### / [SlotPlot](NotUploadedYet) \\
Edited Seurat::DotPlot to split the plot into groups of genes for better visual clarity.

### / [LimPlot](R/LimPlot.R) \\
LimPlot fixes Seurat's DimPlot, which by default requires you specify label = T, label.box = T, repel = T and +NoLegend(). LimPlot has all of these by default, saving literally 10s of seconds.

### / [protanCols](R/protanCols.R) \\
Colourblind-friendly colours for plotting. Currently only 12 but more to be added.

### / [synPanel](R/synPanel.R) \\
Helpful lists of common genes for clustering synovial scRNA-seq data

### / [FindNeighbours](R/FindNeighbours.R) \\
Fixes Seurat's "FindNeighbors", which is spelled incorrectly (masks Seurat::FindNeighbors but with British English Spelling)
