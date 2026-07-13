# Fetch orthologs from InParanoiDb

This function retrieves the orthologs from [InParanoiDB
9](https://inparanoidb.sbc.su.se/) database for 640 species given set of
Ensembl Gene IDs or uniprot IDs usinf database API.

## Usage

``` r
searchIpDb(geneID, ..., idtype = "ensembl")
```

## Arguments

- geneID:

  Gene ID of *Plasmodium falciparum* or VEupathDB enlisted organisms
  that are also covered by InParanoiDB. If providing uniprot ID, set
  idtype="uniprot" to prevent ID conversion.

- ...:

  Additional arguments that can be passed to
  [`toGeneid()`](https://rohit-satyam.github.io/plasmoRUtils/reference/toGeneid.md)
  function. This comes in very handy if you are working with parasite
  gene IDs other than *Plasmodium*.

- idtype:

  Set this to "uniprot" if using uniprot IDs directly and GenID to
  uniprot ID conversion is not required.

## Value

A data frame, containing 10 columns.

## Details

To view list of species covered by InParanoiDB 9, use
[`listipdb()`](https://rohit-satyam.github.io/plasmoRUtils/reference/listipdb.md)
function.

## See also

[`listipdb`](https://rohit-satyam.github.io/plasmoRUtils/reference/listipdb.md)

## Examples

``` r
if (FALSE) { # \dontrun{
df <-  searchIpDb( c("PF3D7_0807800", "PF3D7_1023900"))
df <- searchIpDb( c("C5LD32", "A5KAC7"),idtype = "uniprot")
} # }
```
