## luktils (luke-tils, like Luke Utils) is a package made to make my life easier.

<sub><sup>Yes that's how you say it. It seemed like a good idea at the time.</sub></sup>

There's no guarantee it'll make your life easier.

<sub><sup>Or mine, for that matter.</sub></sup>

scRNA-seq analysis, plotting, and dataviz.

<sub><sup>And it doesn't even do that particularly well.</sub></sup>


## Current Functions:
### / subToMain \\
This is a quick function to move the active idents from an object subclustered new labels to a larger object by matching cell names.

### / avgHeatmap \\
This is a heatmap plotting that uses pheatmap and AggregateExpression to make averaged expression heatmaps per cluster, instead of the default DoHeatmap

### / freqCharts \\
This makes 3 charts which show one metadata against another, both in raw numbers and in frequency (scaled from 0-1) to show, for example, abundance of clusters per experimental group.

### / lUMAP \\
lUMAP (El-UMAP or Loo-MAP, as in "Luke UMAP") uses a previously generated UMAP graph but changes the plot size, point circling, larger legend, and plot boxed in.

### / protanCols \\
Colourblind-friendly colours for plotting. Currently only 12 but more to be added.
