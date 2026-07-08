# Get effective gene lengths form TPM calculation

This function provides effective length (sum of lengths of exons) of the
genes for calculating TPM values.

## Usage

``` r
getEffLen(gtf = NULL, format = "gff3")
```

## Arguments

- gtf:

  Provide path or URL to GTF file.

- format:

  Format of the feature file i.e. "gtf" or "gff3". Default "gff3".

## Value

df numeric. This function returns a data frame with 2 columns: "GeneID",
"Length".

## Examples

``` r
if (FALSE) { # \dontrun{
baseurl <- "https://plasmodb.org/common/downloads/release-68/"
getEffLen(paste0(baseurl,
"Pfalciparum3D7/gff/data/PlasmoDB-68_Pfalciparum3D7.gff"))

OR

getEffLen("/data/PlasmoDB-67_Pfalciparum3D7.gtf")
} # }
```
