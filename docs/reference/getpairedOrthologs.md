# Fetch paired orthologs

This function provides ability to query InParanoiDB 9 and OrthoMCL 7 to
get paired orthologs between two species of interest.

## Usage

``` r
getpairedOrthologs(
  from,
  to,
  db = c("ipdb", "orthomcl"),
  customFields = NULL,
  transform = TRUE,
  api_key = Sys.getenv("VEUPATHDB_API_KEY")
)
```

## Arguments

- from:

  ID of query organism for which orthologs are to be fetched. To view
  organisms indexed in the database and their respective IDs use
  [`listipdb()`](https://rohit-satyam.github.io/plasmoRUtils/reference/listipdb.md)
  or
  [`listOrthomcl()`](https://rohit-satyam.github.io/plasmoRUtils/reference/listOrthomcl.md).

- to:

  ID of target organism against which will be queried for orthologs. To
  view organisms indexed in the database and their respective IDs use
  [`listipdb()`](https://rohit-satyam.github.io/plasmoRUtils/reference/listipdb.md)
  or
  [`listOrthomcl()`](https://rohit-satyam.github.io/plasmoRUtils/reference/listOrthomcl.md).

- db:

  Define database to be queried. Possible values: "orthomcl","ipdb".

- customFields:

  Additional field to be fetched from OrthoMCL 7. "primary_key" and
  "target_id" are mandatory fields. Additional popular fields include:
  "group_name","product","source_id","num_core","num_peripheral","length","sequence","taxon_name","abbreviation","core_peripheral","ec_numbers","pfam_domains"
  etc.For more fields refer to OrthoMCL REST query builder.

- transform:

  Logical. In case of InParanoiDB, transform collapses the orthologs by
  Group IDs so that you have unique rows. For OrthoMCL, it performs the
  opposite function and separates the rows so that you have only one ID
  in each column. When set to FALSE, the function returns the ortholog
  query results in the raw form provided by the database.

- api_key:

  Character. VEuPathDB API key used for authentication. By default, the
  value is obtained from the `VEUPATHDB_API_KEY` environment variable.

## Value

df This function returns a dataframe of orthologs and additional
requested fields.

## Examples

``` r
if (FALSE) { # \dontrun{
ids <- listOrthomcl()
df <- getpairedOrthologs(
  from = ids %>% dplyr::filter(Organism == "Plasmodium falciparum 3D7") %>% dplyr::pull(ID),
  to = ids %>% dplyr::filter(Organism == "Toxoplasma gondii ME49") %>% dplyr::pull(ID),
  db = "orthomcl",
  transform = TRUE
)
df <- df %>%
  dplyr::filter(`Target ID` != "No Ortholog Found") ## remove IDs without orthologs
} # }
```
