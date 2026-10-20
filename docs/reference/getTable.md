# Get tables with custom fields

A convenience function to quickly fetch table of Gene IDs, Protein IDs,
Gene Symbols, Annotations and many more columns from database of your
choice such as PlasmoDB, ToxoDB, PiroplasmaDB among other VEuPathDB
pathogen databases.

## Usage

``` r
getTable(
  org,
  db = "toxodb",
  customFields = NULL,
  api_key = Sys.getenv("VEUPATHDB_API_KEY")
)
```

## Arguments

- org:

  Full name of organism of interest as specified in VEuPathDB. To find
  the exact name of the organism, use `listVeupathdb` function.

- db:

  Character Name of the database in which the organism is present. These
  can be one of the following:
  "toxodb","plasmodb","hostdb","amoebadb","cryptodb","fungidb","giardiadb","microsporidiadb","piroplasmadb","trichdb","tritrypdb".

- customFields:

  A vector of custom fields desired to be fetched. "primary_key" is
  mandatory field. Use
  [`getTableAttributes()`](https://rohit-satyam.github.io/plasmoRUtils/reference/getTableAttributes.md)
  to browse the attributes available for the selected database and
  record type. If `NULL`, a default set of attributes is returned.

- api_key:

  VEuPathDB API key.

## Value

A data frame, containing "Gene ID", "Product Description", "Gene
Strand", "Gene Name or Symbol", "Previous ID(s)", "Entrez Gene ID",
"UniProt ID(s)", "Protein Length", "TM Domains" and "SignalP Peptide"
for all the genes present in the organism of interest.

## Details

With release 71 of VEuPathDB, only users with subscription and API key
can access the database programmatically. Use
`usethis::edit_r_environ(scope = "user")` to add VEuPathDB API using
`VEUPATHDB_API_KEY` variable.

## See also

[`getTableAttributes()`](https://rohit-satyam.github.io/plasmoRUtils/reference/getTableAttributes.md)
for browsing the attributes available for a given record type.

## Examples

``` r
if (FALSE) { # \dontrun{

df <- getTable(
  org = "Plasmodium falciparum 3D7",
  db = "plasmodb", api_key = Sys.getenv("VEUPATHDB_API_KEY")
)

df <- getTable(
  org = "Plasmodium falciparum 3D7",
  db = "plasmodb", api_key = Sys.getenv("VEUPATHDB_API_KEY"),
  customFields = c("primary_key", "protein_sequence")
)
} # }
```
