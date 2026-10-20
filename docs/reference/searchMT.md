# Fetch data tables from malaria.tools database

This function retrieves data from malaria.tools and generates a
dataframe containing the Stage of Parasite in which the gene is highly
expressed.

## Usage

``` r
searchMT(geneID)
```

## Arguments

- geneID:

  A character vector of Gene IDs of Plasmodium falciparum or Plasmodium
  berghi.

## Value

A plot (or data thereof) of TPM values across multiple stages of
parasite.

## Examples

``` r
if (FALSE) { # \dontrun{
  geneID <- c("PBANKA_0100600", "PBANKA_0102900", "PF3D7_0102900")
  ## To get condition specificity and tissue specificity data
  res <- searchMT(geneID = geneID)
} # }
```
