# Get effective gene lengths form TPM calculation

This function provides effective length (sum of lengths of exons) of the
genes for calculating TPM values. Calculates the union of all annotated
exons across every transcript associated with a gene.

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

getEffLen("/data/PlasmoDB-67_Pfalciparum3D7.gtf")
} # }
```
