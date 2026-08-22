# Fetch articles from PubMed database

This function searches the PubMed corpus for the articles that contains
your Gene ID of interest.

## Usage

``` r
searchPM(
  geneID,
  org = "Plasmodium falciparum",
  query = NULL,
  from = 2010,
  to = 2026,
  verbose = TRUE
)
```

## Arguments

- geneID:

  Character vector of Gene IDs.

- org:

  Scientific name of the organism. Default *Plasmodium falciparum*

- query:

  String of user defined custom queries. If you wish to pass your own
  combination od terms use this argument.

- from:

  To define the start year for querying the articles.

- to:

  To define the end year for querying the articles.

- verbose:

  Disable to turn off the messages printed by the function.

## Value

A data frame, containing 9 columns: "pmid" "doi" "title" "year" "month"
"day" "jabbrv" "journal" "GeneID" .

## Examples

``` r
if (FALSE) { # \dontrun{
test <- searchPM(geneID = c("PF3D7_0420300","PF3D7_0621000"))
} # }
```
