# Quickly fetch data from MPMP database

This function provides ability to fetch pathways data from MPMP. This
data can then be modified and used for MPMP pathway enrichment analysis.

## Usage

``` r
getMpmp(url)
```

## Arguments

- url:

  URL of the pathway of interest.

## Value

df This function returns a dataframe of containing the gene ID and
Annotations fetched from MPMP database.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- getMpmp("http://mpmp.huji.ac.il/maps/HNE_prot.html")
df <- getMpmp("http://mpmp.huji.ac.il/maps/14-3-3prot.html")
} # }
```
