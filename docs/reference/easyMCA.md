# Download data from the Malaria Cell Atlas

Download an MCA dataset archive, or read a locally available ZIP
archive, and import the requested metadata or expression table.

## Usage

``` r
easyMCA(url, type = "data", timeout = 300)
```

## Arguments

- url:

  URL of an MCA ZIP archive (for example, one returned by
  [`listMCA()`](https://rohit-satyam.github.io/plasmoRUtils/reference/listMCA.md)),
  or the path to a local ZIP archive.

- type:

  Character string specifying the table to import: `"data"` (cell
  metadata; default), `"exp"` (processed expression values), or `"raw"`
  (raw counts).

- timeout:

  Download timeout in seconds (default: 300). This applies only to
  remote URLs; the previous R timeout option is restored on exit.

## Value

A `data.frame` containing the selected table, with the first CSV column
used as row names. For expression tables, genes are expected to occupy
rows and cells columns, as provided by the MCA download.

## Details

Only the requested CSV file is extracted, into an isolated temporary
directory. Temporary files are removed when the function exits,
including if an error occurs. Local input ZIP archives are never
deleted. The supported archive members end in `-data.csv`, `-exp.csv`,
or `-raw.csv` (also accepting underscore separators and `.csv.gz`).

## Examples

``` r
if (FALSE) { # \dontrun{
url <- "https://www.malariacellatlas.org/downloads/pf-ch10x-set4-biorxiv.zip"

metadata <- easyMCA(url, type = "data")
expression <- easyMCA(url, type = "exp")
raw_counts <- easyMCA(url, type = "raw")

# An already downloaded archive can also be used offline.
metadata <- easyMCA("path/to/pf-ch10x-set4-biorxiv.zip", "data")

# When creating a Seurat object, check that cells match and
# convert the counts data frame to a numeric matrix.
stopifnot(all(colnames(raw_counts) %in% rownames(metadata)))
seurat_obj <- Seurat::CreateSeuratObject(
  counts = as.matrix(raw_counts),
  meta.data = metadata[colnames(raw_counts), , drop = FALSE],
  project = "MCA"
)
} # }
```
