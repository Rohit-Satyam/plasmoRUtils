# plasmoRUtils

A convenience function to quickly convert the Old Pf IDs, Uniprot or
Entrez Ids to Ensembl IDs using PlasmoDB Release 68 Annotation data.It
also provides description and gene symbol for input Ids.

## Usage

``` r
toPfid(inputid, from = "", to = "")
```

## Arguments

- inputid:

  A character vector of IDs. Can be Ensembl, Uniprot, Entrez or or old
  Pf ids.

- from:

  To describe the type of Input ID. Possible values: "old", "uniprot".
  "entrez", "ensembl"

- to:

  To describle the type of output ID desired. Possible values:
  "emsembl".

## Value

A data frame, containing PFIDs, gene description and gene Symbols.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- toPfid(c("PF3D7_0420300", "PF3D7_0621000"), from="ensembl")
} # }
```
