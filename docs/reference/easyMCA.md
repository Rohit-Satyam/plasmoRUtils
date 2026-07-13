# Fetch data from Malaria Cell Atlas

This function fetches and loads the metadata and expression matrices
from the desired data sets available on Malaria Cell Atlas (MCA)
Database.

## Usage

``` r
easyMCA(url, type = "data")
```

## Arguments

- url:

  A url of dataset from `listMCA` function.

- type:

  Type of data to be fetched. Use "exp" to fetch normalized and scaled
  values, use "raw" to get raw counts and use "data" to get the metadata
  of the dataset.

## Value

df A dataframe of requested data type from MCA.

## Examples

``` r
if (FALSE) { # \dontrun{
  url <- "https://www.malariacellatlas.org/downloads/pf-ch10x-set4-biorxiv.zip"
  # Use the function to read metadata, expression, or raw data
  metadata <- easyMCA(url, type = "data")
  expression <- easyMCA(url, type = "exp")
  raw_counts <- easyMCA(url, type = "raw")

  ## make Seurat Object easily now

  testmca <- Seurat::CreateSeuratObject(counts = raw_counts, meta.data = metadata, project = "MCA")
} # }
```
