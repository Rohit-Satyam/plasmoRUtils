# Calculate TPM values from count data

This function provides ability to compute quick TPM values.

## Usage

``` r
easyTPM(counts, featureLength)
```

## Arguments

- counts:

  Count matrix containing raw counts.

- featureLength:

  Effective length of genes generated from `getEfflen`.

## Value

df numeric. This function returns a dataframe of TPM normalized counts
and final column containing feature length.

## Examples

``` r
if (FALSE) { # \dontrun{
## Effective length of the gene
gene_info <- data.frame(GeneID = c("gene1", "gene2", "gene3"), Length = c(1000, 1500, 2000))
count_matrix <- matrix(c(10, 20, 30, 40, 50, 60),
  nrow = 3, ncol = 2,
  dimnames = list(c("gene1", "gene2", "gene3"), c("sample1", "sample2"))
)
test <- easyTPM(count_matrix, gene_info)
} # }
```
