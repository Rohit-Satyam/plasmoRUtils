# Pie chart to summarize searchMT results

This function make a donut chart to represent distribution of input gene
IDs across different stages of Malaria given a result object from
`plotTissueSpecific` function.

## Usage

``` r
easyPie(df, col = "Tissue Specificity")
```

## Arguments

- df:

  A dataframe obtained from plotTissueSpecific(returnData=TRUE).

- col:

  Column to plot as donut chart. Default: "Tissue Specificity"

## Value

A plot (or data thereof) of domains present in the list of gene IDs.

## Examples

``` r
if (FALSE) { # \dontrun{
  geneID <- c("PBANKA_0100600", "PBANKA_0102900", "PF3D7_0102900")
  ## To get Plot similar to malaria.tools
  res <- searchMT(geneID = geneID)
  res %>% easyPie()
} # }
```
