# Get pre-configured tables from OrthoMCL-DB

A convenience function to quickly fetch preconfigured table of from
OrthoMCL-DB. This function can be used to map old OG IDs to new and
vice-versa given a vector of OG IDs. If no IDs are provided, entire data
tables fetched from the database will be returned, which can be
manipulated or processed further as per user need.

## Usage

``` r
getPreconfiguredTableOrthomcl(
  og_ids = NULL,
  customField = "previousGroups",
  timeout = 60,
  max_tries = 3
)
```

## Arguments

- og_ids:

  Orthogroup IDs.

- customField:

  Preconfigured tables to be downloaded from OrthoMCL.

- timeout:

  Timeout time in seconds

- max_tries:

  Maximum retries allowed.

  - *"TaxonCounts"*: Phlyetic Distribution of Proteins

  - *"PhyleticDistributionCounts"*: Phyletic Distribution clade
    abundance

  - *"previousGroups"*: All previous groups."Default"

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
df <- getPreconfiguredTableOrthomcl(customField="SimilarGroup")
} # }
```
