# Calculate transcripts per million (TPM) from counts

Normalise a raw count matrix using matched gene lengths. Gene lengths
are aligned by gene identifier before calculating reads per kilobase and
scaling each sample to one million TPM.

## Usage

``` r
easyTPM(counts, featureLength)
```

## Arguments

- counts:

  Count matrix containing raw counts.

- featureLength:

  Effective length of genes generated from
  [`getEffLen()`](https://rohit-satyam.github.io/plasmoRUtils/reference/getEffLen.md).

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
  dimnames = list(c("gene2", "gene1", "gene3"), c("sample1", "sample2"))
)
test <- easyTPM(count_matrix, gene_info)
} # }
```
