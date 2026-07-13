# Making plots similar to malaria.tools

This function retrieves data from malaria.tools and generates expression
value plots (in TPM) similar to those produced by the website. Use this
function to create publication-ready plots.

## Usage

``` r
plotTissueSpecific(geneID, returnData = FALSE, plotify = FALSE)
```

## Arguments

- geneID:

  Gene ID of Plasmodium falciparum or Plasmodium berghi.

- returnData:

  Logical. Use true to return dataframe used for making plots.

- plotify:

  To make plots interactive using plotly.

## Value

A plot (or data theirof) of TPM values across multiple stages of
parasite.

## Examples

``` r
if (FALSE) { # \dontrun{
  geneID <- c("PBANKA_0100600", "PBANKA_0102900", "PF3D7_0102900")
  ## To get Plot similar to malaria.tools
  res <- plotTissueSpecific(geneID = "PBANKA_0100600")
} # }
```
