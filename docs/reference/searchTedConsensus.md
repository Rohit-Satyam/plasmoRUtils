# Fetch protein domains from TED database

A convenience function to quickly access The Encyclopedia of Domains
(TED) database and fetch domain boundary information for given Uniprot
IDs. For information about the column names the users are requested to
refer to the TED database at https://ted.cathdb.info/.

## Usage

``` r
searchTedConsensus(
  uniprotid = "",
  returnCATHdesc = TRUE,
  timeout = 60,
  max_tries = 3
)
```

## Arguments

- uniprotid:

  A character vector of uniprot IDs.

- returnCATHdesc:

  Logical. Set this on to get the description of the CATH ID from CATH
  database.

- timeout:

  Timeout time in seconds

- max_tries:

  Maximum retries allowed.

## Value

A data frame, domain boundaries and other information provided by TED.
For details visit the TED database.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- searchTedConsensus(
c("Q7K6A1","Q8IAP8","C0H4D0","C6KT90","Q8IBJ7"),
returnCATHdesc=FALSE)
} # }
```
