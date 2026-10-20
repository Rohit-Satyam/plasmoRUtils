# Plot GO enrichment results

Plots the most significant GO terms from
[`easytopGO()`](https://rohit-satyam.github.io/plasmoRUtils/reference/easytopGO.md)
or a VEuPathDB enrichment table. Input p-values are sorted before
selecting the top terms. Redundant rows are identified using GO.ID when
available, otherwise the term description; this is a display-only
reduction, not a change to enrichment.

## Usage

``` r
easyGOPlot(
  res,
  title = "GO Biological processes",
  limit = 20L,
  desc = "Term",
  genecounts = "Significant",
  sortby = NULL,
  min_genes = NULL
)
```

## Arguments

- res:

  A data.frame of GO enrichment results.

- title:

  Plot title.

- limit:

  Maximum number of terms to display (positive integer).

- desc:

  Name of the term-description column (e.g. "Term" or "Name").

- genecounts:

  Name of the gene-count column (e.g. "Significant" or "Result count").
  For KS, Significant may refer to a nominal gene cutoff.

- sortby:

  Name of a numeric GO p-value column. When NULL (default), uses the
  first available among padj, fisher, ks, P-value and p-value.

- min_genes:

  Optional minimum count for plotting (NULL means no filter). Applied
  only to the displayed terms, not before GO enrichment/FDR testing.

## Value

A ggplot2 object. Empty results produce a labelled empty plot.

## Examples

``` r
if (FALSE) { # \dontrun{
gores <- easytopGO(
  geneID = down_padj, bkggset = background,
  gaf = "annotation.gaf.gz", useBiomart = FALSE,
  useGAF = TRUE, stats = "fisher"
)
easyGOPlot(gores, sortby = "fisher", min_genes = 3)
# For KS output:
# easyGOPlot(ks_result, sortby = "ks")
# For a VEuPathDB result with numeric p-values:
# easyGOPlot(veupath, desc = "Name", genecounts = "Result count",
#            sortby = "P-value")
} # }
```
