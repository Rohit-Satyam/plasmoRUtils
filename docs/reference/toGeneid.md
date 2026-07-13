# Convert Other IDs to Ensembl gene IDs

A convenience function to quickly convert the Uniprot or Entrez Ids to
Ensembl gene IDs, using VEuPathDB specialized databases.It also provides
description and gene symbol for input Ids.

## Usage

``` r
toGeneid(
  inputid,
  from = "",
  to = "",
  org = "Plasmodium falciparum 3D7",
  db = "plasmodb",
  ...
)
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

- org:

  Organism to which the IDs belongs to.Possible values: toxodb,
  plasmodb, hostdb, amoebadb, cryptodb, fungidb, giardiadb,
  microsporidiadb, piroplasmadb, trichdb, tritrypdb.

- db:

  Database in which organism is present.

- ...:

  Additional arguments that can be passed to the `getTable` function.

## Value

A data frame, containing Gene IDs, gene description and gene Symbols and
more.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- toGeneid(
c("PF3D7_0420300", "PF3D7_0621000"),
      from="ensembl")
} # }
```
