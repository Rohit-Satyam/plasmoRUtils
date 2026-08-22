# Fetch articles from Research Square database

This function searches the Research Square corpus for the articles that
contains your Gene ID/ Gene name of interest.

## Usage

``` r
searchRS(
  geneID,
  org,
  gene_aliases = character(),
  org_aliases = character(),
  max_pages = 5,
  sleep_sec = 0.2,
  status = "all",
  article_type = "Research Article"
)
```

## Arguments

- geneID:

  Gene ID such as PF3D7_XXXX or Gene Symbol.

- org:

  Scientific name of the organism. Eg. *Plasmodium falciparum*

- gene_aliases:

  Gene aliases, short form or alternate names.

- org_aliases:

  Common name or aliases to filter irrelevant articles.

- max_pages:

  Maximum number of pages to parse.

- sleep_sec:

  Sleep time between subsequent queries.

- article_type:

  The type of Article to return. Possible values: "Research Article",
  "Systematic Review", "Method Article", "Short Report", "Case Report",
  "Data Note", "Video".

## Value

A data frame, containing 10 columns such as Research square ID, Title,
Authors, Date, Status, Title of the journal, Article type, URLs, Gene ID
present (logical), Organism name present (logical).

## Examples

``` r
if (FALSE) { # \dontrun{

hits <- searchRS(
geneID = "PfMEDA1",
org = "Plasmodium falciparum",
max_pages = 10
)

hits <- searchRS(
geneID = "PfAP2-I",
org = "Plasmodium falciparum",
gene_aliases = c(
  "Pf AP2-I",
  "pfap2-i",
  "PF3D7_1007700",
  "Apetala 2 Invasion",
  "AP2-I transcription factor"
),
org_aliases = c(
  "P. falciparum",
  "3D7",
  "malaria parasite"
),
max_pages = 10
)
} # }
```
