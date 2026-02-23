# annotated-volcano-plot

A lightweight R function to create publication-ready volcano plots from differential expression results (e.g., DESeq2, edgeR, limma). The function supports customizable cutoffs, gene labeling with ggrepel, legend counts, and intelligent x-axis clipping with visual marking of out-of-range points.

## Features

* Publication-ready volcano plots using ggplot2
* Automatic classification of genes (Up-regulated / Down-regulated / Non-significant)
* Custom log2FC and adjusted p-value cutoffs
* Optional gene labeling with ggrepel
* Intelligent x-axis clipping (points outside limits are clamped and visually marked)
* Legend with category counts (N)
* Flexible color customization
* Automatic handling of NA and zero adjusted p-values
* Optional file export via ggsave

## Requirements

This function requires the following R packages:

* ggplot2
* ggrepel

Install them with:

```r
install.packages(c("ggplot2", "ggrepel"))
```

## Installation

Clone the repository:

```bash
git clone https://github.com/michelemarziliano/annotated-volcano-plot.git
cd annotated-volcano-plot
```

Then source the function in R:

```r
source("annotated_volcano_plot.R")
```

## Input Data Format

The function expects a data.frame containing differential expression results with:

* A column named `log2FoldChange`
* A column named `padj` (adjusted p-values)
* Row names corresponding to gene IDs

Example structure:

```r
head(dea_res)
#              log2FoldChange      padj
# GENE1              1.25        0.001
# GENE2             -0.87        0.12
# GENE3              2.10        0.0005
```

## Quick Start Example

```r
library(ggplot2)
library(ggrepel)

source("annotated_volcano_plot.R")

thr_lfc <- 1
thr_padj <- 0.05

annotated_volcano_plot(
  dea_res = volcanodata,
  main = paste0(
    "Volcano plot of NR vs R\n",
    "with ABS(Log2FC) > ", thr_lfc,
    " and adjusted P-value < ", thr_padj
  ),
  cutoff_log2FC = thr_lfc,
  cutoff_padj = thr_padj,
  genes_to_label = genes_to_label,
  col_up_genes = "red",
  col_down_genes = "blue",
  col_other_genes = "grey40",
  file_path = "volcano_NR_vs_R.png",
  xlim_range = c(
    -max(abs(volcanodata$log2FoldChange)),
     max(abs(volcanodata$log2FoldChange))
  ),
  width_in = 6,
  height_in = 6,
  show_legend = TRUE
)
```

## Function Reference

### annotated_volcano_plot()

```r
annotated_volcano_plot(
  dea_res,
  main = "",
  cutoff_log2FC = 1,
  cutoff_padj = 0.05,
  deg_list = NULL,
  genes_to_label = NULL,
  col_up_genes,
  col_down_genes,
  col_other_genes,
  file_path = NULL,
  report_cutoffs = FALSE,
  descr = "",
  xlim_range = NULL,
  width_in = 8,
  height_in = 8,
  show_legend = TRUE
)
```

### Arguments

**dea_res**
Data frame containing differential expression results. Must include the columns `log2FoldChange` and `padj`. Row names must be gene IDs.

**main**
Main plot title (character string).

**cutoff_log2FC**
Numeric threshold for absolute log2 fold change (|log2FC|).

**cutoff_padj**
Numeric threshold for adjusted p-value.

**deg_list** (optional)
Character vector of gene IDs or logical/integer indexing specifying differentially expressed genes directly, instead of using cutoffs. Genes are classified as Up- or Down-regulated based on the sign of log2FoldChange.

**genes_to_label** (optional)
Character vector of gene IDs to label on the plot. Missing genes are ignored with a warning.

**col_up_genes**
Color for up-regulated genes (any valid R color).

**col_down_genes**
Color for down-regulated genes (any valid R color).

**col_other_genes**
Color for non-significant genes (any valid R color).

**file_path** (optional)
File path (with extension) to save the plot using ggsave. If NULL, the plot is printed instead of saved.

**report_cutoffs**
Logical. If TRUE, cutoff values are automatically appended to the plot title.

**descr**
Optional short description appended to the title (ignored if report_cutoffs = TRUE).

**xlim_range** (optional)
Numeric vector of length 2 specifying x-axis limits. Points outside the range are clipped to the boundary and visually marked with a different shape.

**width_in, height_in**
Plot dimensions in inches when saving the file.

**show_legend**
Logical. If TRUE, the legend is displayed at the bottom; if FALSE, it is hidden.

## Notes

* Adjusted p-values (padj) equal to 0 are internally replaced with `.Machine$double.xmin` to avoid infinite -log10 values.
* NA padj values are automatically excluded from plotting.
* The legend order is fixed (Down-regulated, Non-significant, Up-regulated) even if some categories are absent.
* Clipped points (outside x-axis limits) are visually distinguished using a different shape.

## Output

* Returns: ggplot object (printed if file_path is NULL)
* Optionally saves the plot to disk if file_path is provided

## Example Output

An example plot can be found in the `figures/` directory of this repository.

## License

This project is released under the MIT License.

## Citation

If you use this function in a publication, please cite the GitHub repository:

Marziliano, M. annotated-volcano-plot. GitHub repository.
