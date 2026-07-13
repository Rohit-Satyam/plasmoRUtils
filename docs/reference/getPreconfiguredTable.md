# Get pre-configured tables

A convenience function to quickly fetch preconfigured table of Signal
Peptide ranges, Pathways, Pubmed entries related to genes, Annotations
and etc from database of your choice such as PlasmoDB, ToxoDB,
PiroplasmaDB among other VEuPathDB pathogen databases.

## Usage

``` r
getPreconfiguredTable(org, db = "plasmodb", customField = "Y2hInteractions")
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
  "MassSpecMod", "Epitopes" etc.

## Value

A data frame.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- getPreconfiguredTable(org = "Plasmodium falciparum 3D7",
     db = "plasmodb",customField = "Y2hInteractions")
} # }
```
