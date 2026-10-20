# Performing quick ORA analysis using Fisher or KS and an explicit background

Runs topGO enrichment using an explicit eligible-gene universe. Fisher
tests over-representation of selected genes; KS uses genuine gene-level
p-value ranks for the complete background. Results can be plotted using
[`easyGOPlot()`](https://rohit-satyam.github.io/plasmoRUtils/reference/easyGOPlot.md).

## Usage

``` r
easytopGO(
  geneID,
  bkggset = NULL,
  gaf = "",
  useBiomart = TRUE,
  useGAF = FALSE,
  mart = "protists_mart",
  gset = "pfalciparum_eg_gene",
  algo = "weight01",
  stats = "fisher",
  category = "BP",
  fdr = FALSE,
  correction = "BY",
  cutoff = 0.05,
  go_cutoff = 0.05,
  return_GOdata = FALSE,
  OrgDb = NULL,
  orgdb_keytype = "GID"
)
```

## Arguments

- geneID:

  Uniquely named numeric vector of p-values in, named by gene ID. For
  Fisher, either selected genes or the complete tested vector may be
  supplied; genes with p-value \< `cutoff` are selected. For KS, supply
  the actual p-value for EVERY gene in `bkggset` (not just DEGs).

- bkggset:

  Nonempty character vector of all genes eligible for the analysis,
  typically genes eligible under your RNA-seq selection procedure.
  Required for both tests.

- gaf:

  Path or URL to a GO Annotation File (GAF; .gz is supported). GAF
  column 2 must contain gene IDs corresponding to `bkggset`. When using
  this argument, set useBiomart=FALSE and useGAF=TRUE.

- useBiomart:

  Logical To enable usage of BiomaRt to fetch GO terms. Default: TRUE.

- useGAF:

  Logical To enable usage of custom .gaf file to fetch GO terms.
  Default: FALSE

- mart, gset:

  BioMart name and dataset, respectively.

- algo:

  Argument to specify which algorithm to be used for enrichment by
  topGO. For possible options use
  [`topGO::whichAlgorithms()`](https://rdrr.io/pkg/topGO/man/getSigGroups.html)

- stats:

  Argument to specify the statistical test to be used for enrichment by
  topGO. For possible values use
  [`topGO::whichTests()`](https://rdrr.io/pkg/topGO/man/getSigGroups.html).
  Default: "fisher".

- category:

  Specify the category of Over-representation analysis such as "BP" for
  Biological Process, "MF" for Molecular Function and "CC" for Cellular
  Component Enrichment.

- fdr:

  logical. Perform multiple testing correction testing. Default (FALSE)

- correction:

  Adjustment method accepted by
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Possible
  values: ""holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr".
  Read section 6.2 in topGO documentation before performing correction.
  Correction when using elim and weight is usually not recommended.

- cutoff:

  P-value threshold defining selected genes, default 0.05. For KS this
  changes selected-gene reporting, not the KS statistic.

- go_cutoff:

  Threshold for reporting raw GO p-values or adjusted GO p-values when
  `fdr = TRUE`.

- return_GOdata:

  If TRUE, return the underlying `topGOdata` object.

- OrgDb:

  An OrgDb object loaded by the caller, e.g. after running
  [`library(org.Pfalciparum3D7.eg.db)`](https://rdrr.io/r/base/library.html).
  Passing it selects offline GO annotations instead of the default
  BioMart route.

- orgdb_keytype:

  The key type of the IDs in `bkggset` and `geneID`. Default `"GID"` for
  custom packages built with AnnotationForge; check the keytypes
  available in your particular database.

## Value

A data.frame with GO IDs, terms, numeric `fisher` or `ks` results,
optional `padj`, and associated/selected genes. The `gene_universe`
attribute gives available, feasible and selected-feasible counts. When
`return_GOdata = TRUE`, returns a `topGOdata` object. Full descriptions
of GO terms truncated by GenTable are restored when pRoloc is installed.

## Details

The OrgDb route uses GO (direct annotations), leaving ancestor
propagation to topGO. topGO itself depends on AnnotationDbi, so no
additional installation or declared direct Import is needed for this
route.

## Examples

``` r
if (FALSE) { # \dontrun{
# res: DESeq2 results with Geneid, pvalue, padj and log2FoldChange.
ok <- !is.na(res$padj) & !is.na(res$Geneid)
background.gset <- res$Geneid[ok]
down <- ok & !is.na(res$log2FoldChange) &
  res$padj < 0.05 & res$log2FoldChange < 0
geneList <- stats::setNames(res$padj[down], res$Geneid[down])

# Fisher accepts a DEG-only list; the full universe is background.gset.
gores.fisher <- easytopGO(
  geneID = geneList, bkggset = background.gset,
  gaf = "annotation.gaf.gz", useBiomart = FALSE, useGAF = TRUE,
  stats = "fisher", category = "BP"
)

# Load your custom OrgDb package outside the function.
library(org.Pfalciparum3D7.eg.db)
pfdb <- org.Pfalciparum3D7.eg.db
gores.org <- easytopGO(
  geneID = geneList, bkggset = background.gset,
  OrgDb = pfdb, orgdb_keytype = "GID",
  stats = "fisher", category = "BP"
)
# KS
ok.ks <- ok & !is.na(res$padj)
bg.ks <- res$Geneid[ok.ks]
geneList.ks <- stats::setNames(res$padj[ok.ks], bg.ks)
gores.ks <- easytopGO(
  geneID = geneList.ks, bkggset = bg.ks,
  gaf = "annotation.gaf.gz", useBiomart = FALSE, useGAF = TRUE,
  stats = "ks", category = "BP"
)
} # }
```
