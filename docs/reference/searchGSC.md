# Fetch articles for given gene IDs from Google Scholar

This function searches the Google Scholar corpus recursively for the
articles that contains your Gene ID of interest.

## Usage

``` r
searchGSC(
  geneIDs,
  year_start = NULL,
  year_end = NULL,
  max_pages = 2,
  sleep_secs = 10,
  verbose = TRUE,
  translate = NULL
)
```

## Arguments

- geneIDs:

  Character vector of Gene IDs. If you want to use gene symbols, use
  organism name alongside to avoid articles that might have similar
  abbreviated word. eg. use "AP2-P AND Plasmodium".

- year_start:

  Limit the results to starting year of interest.

- year_end:

  Limit the results to end year of interest.

- max_pages:

  Maximum number of pages to scrap.

- verbose:

  Print warnings.

- translate:

  Translate the paper titles to english or desired language. Use two
  letter code.eg: "fr" for french, "en" for english and "es" for
  spanish.

## Value

A data frame, containing 5 columns: GeneID, Title of the article, Year
of Publication, Url and Authors.

## Details

**Warning**: Scraping Google Scholar is against their [Terms of
Service](https://scholar.google.com/intl/en/scholar/about.html). We
advise users to use this function for querying few IDs (not more than
20) per day. Proceeding with this function may result in your IP being
blocked temporarily.

## Examples

``` r
if (FALSE) { # \dontrun{
## We have a fake ID: PF3D7_0420300OR
res <- searchGSC(
geneIDs=c("PF3D7_0420300 OR MAL4P1.192 OR Q8I1N6 OR PFD0985w","PF3D7_0621000","PF3D7_0420300OR"),
translate = "en",
year_start = 2018, 
year_end   = 2021)

test <- searchGSC(geneID = c("AP2-P AND Plasmodium", "AP2-I"))
} # }
```
