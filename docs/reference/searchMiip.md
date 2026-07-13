# Fetch Protein-Protein interactions from MIIP database

This function retrieves Protein-protein interaction data from [MIIP
database](http://www.hpppi.iicb.res.in/pfnet/).

## Usage

``` r
searchMiip(geneID)
```

## Arguments

- geneID:

  A character vector of Gene IDs of *Plasmodium falciparum*.

## Value

A data frame of Protein protein interaction provided by MIIP database.

## Examples

``` r
if (FALSE) { # \dontrun{
 df <- searchMiip(c("PF3D7_0807800","PF3D7_1023900"))
} # }
```
