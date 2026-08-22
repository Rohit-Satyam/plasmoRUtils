# Get tables with custom fields

A convenience function to quickly fetch table of Gene IDs, Protein IDs,
Gene Symbols, Annotations and many more columns from database of your
choice such as PlasmoDB, ToxoDB, PiroplasmaDB among other VEuPathDB
pathogen databases.

## Usage

``` r
getTable(org, db = "toxodb", customFields = NULL)
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
  mandatory field. Other fields can be supplied and can be chosen from
  (but are not limited to): "organism", "gene_location_text",
  "gene_product", "gene_type", "exon_count", "gene_exon_count",
  "gene_transcript_count", "three_prime_utr_length",
  "five_prime_utr_length", "strand", "is_pseudo", "transcript_length",
  "is_deprecated", "gene_name", "gene_source_id", "transcript_product",
  "protein_length", "chromosome", "location_text", "sequence_id",
  "gene_ortholog_number", "gene_orthomcl_name", "gene_paralog_number",
  "cds_length", "molecular_weight", "isoelectric_point", "tm_count",
  "signalp_peptide", "predicted_go_id_component",
  "predicted_go_component", "predicted_go_id_function",
  "predicted_go_function", "predicted_go_id_process",
  "predicted_go_process", "annotated_go_id_component",
  "annotated_go_component", "annotated_go_id_function",
  "annotated_go_function", "annotated_go_id_process",
  "annotated_go_process", "ec_numbers", "ec_numbers_derived"

## Value

A data frame, containing "Gene ID", "Product Description", "Gene
Strand", "Gene Name or Symbol", "Previous ID(s)", "Entrez Gene ID",
"UniProt ID(s)", "Protein Length", "TM Domains" and "SignalP Peptide"
for all the genes present in the organism of interest.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- getTable(org="Plasmodium falciparum 3D7", db="plasmodb")
} # }
```
