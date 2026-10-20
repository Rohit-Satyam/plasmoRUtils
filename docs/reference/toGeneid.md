# Convert identifiers from one type to another

Converts between current gene IDs, previous gene IDs, UniProt IDs and
Entrez IDs using getTableVeupathdb(). All possible combinations are
allowed.

## Usage

``` r
toGeneid(
  inputid,
  from = "",
  to = "ensembl",
  api_key = Sys.getenv("VEUPATHDB_API_KEY")
)
```

## Arguments

- inputid:

  Character vector of identifiers to convert.Can be vector of mixed
  identifiers.

- from:

  Input identifier type. One of "ensembl", "uniprot", "entrez" or "old".

- to:

  Desired output identifier type. One of "ensembl", "uniprot", "entrez"
  or "old". Defaults to "ensembl".

- api_key:

  VEuPathDB API key.

## Value

A data frame containing Input ID, Output ID, Gene ID and Organism.
Multiple mappings are returned as separate rows.

## Details

Supports mixed-organism identifiers without requiring organism or
database information.

All available gene-identifier associations are retained, including
one-to-many and many-to-many mappings.

With release 71 of VEuPathDB, only users with subscription and API key
can access the database programmatically. Use
`usethis::edit_r_environ(scope = "user")` to add VEuPathDB API using
`VEUPATHDB_API_KEY` variable.

## Examples

``` r
if (FALSE) { # \dontrun{

# Gene IDs to UniProt IDs
df <- toGeneid(
  c("PF3D7_1372300", "PBANKA_1355600"),
  from = "ensembl",
  to = "uniprot"
)

# UniProt IDs to gene IDs
df <- toGeneid(
  c("A0A509APC6", "Q75UY1"),
  from = "uniprot",
  to = "ensembl"
)

# Previous IDs to UniProt IDs
df <- toGeneid(
  "PBANKA_135560",
  from = "old",
  to = "uniprot"
)
} # }
```
