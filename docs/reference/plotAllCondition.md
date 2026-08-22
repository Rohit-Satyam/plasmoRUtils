# Making plots similar to malaria.tools

This function retrieves data from malaria.tools and generates expression
value plots (in TPM) similar to those produced by the website. Use this
function to create publication-ready plots.

## Usage

``` r
plotAllCondition(geneID, returnData = FALSE, plotify = FALSE)
```

## Arguments

- geneID:

  Single Gene ID of Plasmodium falciparum or Plasmodium berghi.

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
  #'   ## To get Plot similar to malaria.tools
  res <- plotAllCondition(geneID = "PBANKA_0100600")
} # }
```
