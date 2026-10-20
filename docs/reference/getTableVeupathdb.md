# Retrieve gene annotations using a mixed list of IDs

Retrieves gene annotations across all databases integrated into
VEuPathDB, without requiring organism or database information.

## Usage

``` r
getTableVeupathdb(
  ids,
  customFields = NULL,
  api_key = Sys.getenv("VEUPATHDB_API_KEY")
)
```

## Arguments

- ids:

  Character vector of gene identifiers.

- customFields:

  Optional character vector of attributes. If NULL, default mapping
  attributes are returned. See getTableAttributes() for available
  attributes.

- api_key:

  VEuPathDB API key.

## Value

A data frame containing the input IDs, matched gene IDs, organism
information and requested annotations.

## Details

Uses the VEuPathDB GeneByLocusTag search.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- getTableVeupathdb(
  ids = c("PF3D7_1372300", "PBANKA_1355600","A0A143ZZR8","W7JSP8" )
)
} # }
```
