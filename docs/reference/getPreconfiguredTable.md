# Get pre-configured tables

A convenience function to quickly fetch preconfigured table of Signal
Peptide ranges, Pathways, Pubmed entries related to genes, Annotations
and etc from database of your choice such as PlasmoDB, ToxoDB,
PiroplasmaDB among other VEuPathDB pathogen databases.

## Usage

``` r
getPreconfiguredTable(
  org = NULL,
  db = "plasmodb",
  customField = "Y2hInteractions",
  listtables = FALSE,
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

- customField:

  Preconfigured table that you wish the fetch. Pass only one value at a
  time from the following: "GeneModelDump", "GeneTranscripts", "Alias",
  "GeneLinkouts", "GeneLocation", "PubMed", "OrthologsLite",
  "LowComplexity", "PdbSimilarities", "3dPreds", "AlphaFoldLinkouts",
  "ProteinProperties", "InterPro", "SignalP", "TMHMM", "ECNumbers",
  "ECNumbersInferred", "protein_length", "chromosome", "location_text",
  "sequence_id", "gene_ortholog_number", "gene_orthomcl_name",
  "gene_paralog_number", "MetabolicPathwaysMPMP", "MetabolicPathways",
  "CompoundsMetabolicPathways", "Y2hInteractions", "MassSpecDownload",
  "MassSpecMod", "Epitopes" etc. You can view list of 58 preconfigured
  tables across VEuPathDB databases using
  `getPreconfiguredTable(listtables = TRUE)`

- listtables:

  To get list of preconfigured tables that can be supplied as an
  argument to `customField`.

- api_key:

  Character. VEuPathDB API key used for authentication. By default, the
  value is obtained from the `VEUPATHDB_API_KEY` environment variable.

## Value

A data frame.

## Details

With release 71 of VEuPathDB, only users with subscription and API key
can access the database programmatically. Use
`usethis::edit_r_environ(scope = "user")` to add VEuPathDB API using
`VEUPATHDB_API_KEY` variable.

## Examples

``` r
if (FALSE) { # \dontrun{
getPreconfiguredTable(listtables = TRUE)
df <- getPreconfiguredTable(
  org = "Plasmodium falciparum 3D7",
  db = "plasmodb", customField = "Y2hInteractions"
)
} # }
```
