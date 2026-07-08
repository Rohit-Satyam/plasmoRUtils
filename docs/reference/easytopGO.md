# Performing quick ORA analysis

A convenience function to quickly perform GO Term Enrichment analysis
using TopGO.The results then can be plotted using `easyGOPlot`

## Usage

``` r
easytopGO(
  geneID,
  bkggset = "",
  gaf = "",
  useBiomart = TRUE,
  useGAF = FALSE,
  mart = "protists_mart",
  gset = "pfalciparum_eg_gene",
  algo = "weight01",
  stats = "ks",
  category = "BP",
  fdr = FALSE,
  correction = "BY"
)
```

## Arguments

- geneID:

  A vector of named p-values. The names should be gene IDs.

- bkggset:

  A character vector of gene IDs that should be used as background.

- gaf:

  A URL or path to the .gaf file obtained from PlasmoDB should you
  choose not to use biomaRt. When using this argument, set
  useBiomart=FALSE and useGAF=TRUE.

- useBiomart:

  Logical To enable usage of BiomaRt to fetch GO terms. Default: TRUE.

- useGAF:

  Logical To enable usage of custom .gaf file to fetch GO terms.
  Default: FALSE

- mart:

  Argument to specify the mart for BiomaRt functions. Default:
  "protists_mart".

- gset:

  Argument to specify the geneset to be used by BiomaRt functions.
  Default: "pfalciparum_eg_gene"

- algo:

  Argument to specify which algorithm to be used for enrichment by
  topGO. For possible options use
  [`topGO::whichAlgorithms()`](https://rdrr.io/pkg/topGO/man/getSigGroups.html)

- stats:

  Argument to specify the statistical test to be used for enrichment by
  topGO. For possible values use
  [`topGO::whichTests()`](https://rdrr.io/pkg/topGO/man/getSigGroups.html).
  Default: "ks".

- category:

  Specify the category of Over-representation analysis such as "BP" for
  Biological Process, "MF" for Molecular Function and "CC" for Cellular
  Component Enrichment.

- fdr:

  logical. Perform multiple testing correction testing. Default (FALSE)

- correction:

  Method to be used to calculate adjusted p-value. Possible values:
  ""holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr". Read
  section 6.2 in topGO documentation before performing correction.
  Correction when using elim and weight is usually not recommended.

## Value

A dataframe of enriched terms with GO description and genes filteres by
uncorrected p-values.

## Examples

``` r
if (FALSE) { # \dontrun{
## making gene list from DESEq2
geneList <- subset(res, regulate=="Up") %>% .$padj
names(geneList) <- subset(res, regulate=="Up") %>% .$Geneid

## background genes will be the genes tested for differential expression
background.gset <- res$Geneid
baseurl <- "https://plasmodb.org/common/downloads/Current_Release/"
url<-paste0(baseurl,"Pfalciparum3D7/gaf/PlasmoDB-68_Pfalciparum3D7_GO.gaf.gz")
gores<-easytopGO(geneID = geneList,useGAF = TRUE,useBiomart = FALSE,gaf=url,
bkggset = background.gset, category = "BP", stats = "ks")
} # }
```
