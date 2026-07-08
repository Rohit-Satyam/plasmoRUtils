# Plotting domains on protein body

This function retrieves data from malaria.tools and generates a
dataframe containing the Stage of Parasite in which the gene is highly
expressed.

## Usage

``` r
plotDomain(
  geneID,
  mart = "protists_mart",
  gset = "pfalciparum_eg_gene",
  input = "ensembl_gene_id",
  fetchid = "uniprotsptrembl",
  returnData = FALSE
)
```

## Arguments

- geneID:

  A character vector of Gene IDs of Plasmodium falciparum or Plasmodium
  berghi.Remove version from the gene ids.

- mart:

  Name of Ensembl Biomart. Default: "protists_mart"

- gset:

  Gene-set of organism. This should be changed if you are dealing with
  other Plasmodium species. Default: "pfalciparum_eg_gene"

- input:

  Input id type. The gene IDs from PlasmoDB are "ensembl_gene_id".

- fetchid:

  Desired output ids

- returnData:

  To return converted gene ids with domain information fetched from
  uniprot.

## Value

A plot (or data thereof) of domains present in the list of gene IDs.

## Examples

``` r
if (FALSE) { # \dontrun{
  ## Search proper mart
  ## View(listDatasets(biomaRt::useEnsemblGenomes(biomart = "protists_mart")))
  ## To get domain information from uniprot and prepare publication ready figures
  plot <- plotDomain(geneID = c("PF3D7_0518900", "PF3D7_0602800", "PF3D7_0624600"))

  ## Currently pberghei doesn't work, see issue: https://github.com/grimbough/biomaRt/issues/110
  plot <- plotDomain(geneID = c("PBANKA_0100600", "PBANKA_0102900"), gset = "pberghei_eg_gene")

  ## Change plot domain colors, if desired. Say you have 15 domains in all 3 proteins combined
  palette <- randomcoloR::distinctColorPalette(15)
  plot + scale_fill_manual(values = palette)
} # }
```
