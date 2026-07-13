# Fetch data tables from PhenoPlasm database

This function searches the Phenotypes of the gene IDs in Phenoplasm
database and enables users to fetch sub-tables such as Disruptability
and Mutant phenotypes.

## Usage

``` r
searchPhPl(geneID = "", org = "pf", fetch = 1)
```

## Arguments

- geneID:

  Character vector of Gene IDs.

- org:

  Abbreviation of the organism. Default "pf"

- fetch:

  Numeric. Use 1 to fetch the "Disruptability" table and 2 to fetch
  "Mutant phenotypes" table.

  - pb: *Plasmodium berghii*

  - pk: *Plasmodium knowlesi*

  - pf: *Plasmodium falciparum*

  - pc: *Plasmodium chabaudi*

  - py: *Plasmodium yoelii*

## Value

A data frame.

## Examples

``` r
if (FALSE) { # \dontrun{
## get phenotype for few genes in plasmodium falciparum
df <- searchPhPl(geneID = c("PF3D7_0420300","PF3D7_0621000","PF3D7_0523800"), org="pf")
df <- searchPhPl(geneID = c("PF3D7_0420300","PF3D7_0621000","PF3D7_0523800"), org="pf", fetch=2)

} # }
```
