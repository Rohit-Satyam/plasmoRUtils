# Get pre-configured tables from OrthoMCL-DB

A convenience function to quickly fetch preconfigured tables from
OrthoMCL-DB. This function can be used to map old OG IDs to new and
vice-versa given a vector of OG IDs. If no IDs are provided, the entire
data table fetched from the database is returned and a copy of it is
cached to reduce the next query time.

## Usage

``` r
getPreconfiguredTableOrthomcl(
  og_ids = NULL,
  customField = "previousGroups",
  api_key = Sys.getenv("VEUPATHDB_API_KEY"),
  timeout = 60,
  max_tries = 3,
  refresh = FALSE
)
```

## Arguments

- og_ids:

  Orthogroup IDs.

- customField:

  Preconfigured table to be downloaded from OrthoMCL.

- api_key:

  Character. VEuPathDB API key used for authentication. By default, the
  value is obtained from the `VEUPATHDB_API_KEY` environment variable.

- timeout:

  Timeout time in seconds.

- max_tries:

  Maximum retries allowed.

- refresh:

  Logical. If `TRUE`, ignores any locally cached copy and downloads the
  requested table again from OrthoMCL. Defaults to `FALSE`.

  - *"TaxonCounts"*: Phyletic Distribution of Proteins

  - *"PhyleticDistributionCounts"*: Phyletic Distribution clade
    abundance

  - *"previousGroups"*: All previous groups. "Default"

  - *"GroupStat"*: Group Statistics

  - *"SimilarGroup"*: Similar Groups

  - *"EcNumber"*: Summary of EC Numbers

  - *"PFams"*: Summary of Pfam domains

  - *"ProteinPFams"*: PFam Architecture of Each Protein

  - *"KeywordFrequency"*: Keyword Frequency

## Value

A data frame.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- getPreconfiguredTableOrthomcl("OG3_10277")
df <- getPreconfiguredTableOrthomcl(customField = "SimilarGroup")
} # }
```
