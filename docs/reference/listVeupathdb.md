# List genomes and metadata in VEuPathDB

A convenience function to quickly fetch table of genomes and their
associated metadata from VEupathDB.

## Usage

``` r
listVeupathdb(customFields = NULL)
```

## Arguments

- customFields:

  A vector of custom fields desired to be fetched. "primary_key" is
  mandatory field. Other fields can be supplied and can be chosen from
  (but are not limited to): "annotation_source", "annotation_version",
  "arraygenecount", "chipchipgenecount", "chromosomeCount",
  "codinggenecount", "communitycount", "contigCount", "ecnumbercount",
  "estcount", "genecount", "genecount_number", "genome_source",
  "genome_version", "gocount", "is_in_apollo", "is_reference_strain",
  "megabps", "ncbi_tax_id", "ncbi_taxon_url", "organism",
  "organism_full", "orthologcount", "othergenecount", "popsetcount",
  "project_id", "proteomicscount", "pseudogenecount", "rnaseqcount",
  "rtpcrcount", "snpcount", "species", "species_ncbi_tax_id",
  "species_ncbi_taxon_url", "supercontigCount", "tfbscount",
  "URLcdsFasta", "URLGenomeFasta", "URLgff", "URLproteinFasta",
  "URLtranscriptFasta". For more fields, refer to the VEuPathDB
  Documentation

## Value

A data frame containing information about genomes present in VEuPathDB
and their attributes.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- listVeupathdb()
df <- listVeupathdb(customFields=c("species", "project_id"))
} # }
```
