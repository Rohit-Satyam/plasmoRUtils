# plasmoRUtils

Hidden function called by
[`easytopGO()`](https://rohit-satyam.github.io/plasmoRUtils/reference/easytopGO.md)
to run topGO.

## Usage

``` r
.usetopGO(
  stats = "fisher",
  category = "BP",
  geneID,
  bkggset,
  gene_2_GO = NULL,
  algo = "weight01",
  fdr = FALSE,
  correction = "BY",
  cutoff = 0.05,
  go_cutoff = 0.05,
  return_GOdata = FALSE,
  OrgDb = NULL,
  orgdb_keytype = "GID"
)
```
