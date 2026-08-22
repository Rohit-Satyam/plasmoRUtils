# Fetch Protein-protein interaction for given gene IDs from Hit Predict database

This function searches the Hitpredict database to retrieve the
Experimental Protein-Protein Interaction data.

## Usage

``` r
searchHP(
  geneID,
  taxid = "36329",
  uniprotToGID = TRUE,
  timeout = 60,
  max_tries = 3
)
```

## Arguments

- geneID:

  Single gene ID.

- taxid:

  Taxon ID of the organism of interest. Default: 36329. If taxon id of
  organism is not known set this to NULL.

- uniprotToGID:

  To convert Uniprot ID to gene ID. Set TRUE for *Plasmodium* geneIDs
  only.

- timeout:

  Seconds after which timeout must be considered.

- max_tries:

  Maximum retries to access the database.

## Value

A data frame, containing 11 columns: "Interaction", "Interactor",
"Name", "Experiments", "Category", "Method.Score", "Annotation.Score",
"Interaction.Score", "Confidence", "QueryID", "ensembl_gene_id".

## Examples

``` r
if (FALSE) { # \dontrun{
test <- searchHP("Q8I1Q4")

## To use it for other organism, turn off uniprotToGID and provide taxid of the organism
test <- searchHP("BRCA1",taxid = "9606" , uniprotToGID = FALSE)
} # }
```
