# Fetch Minor intron tables from MiDB

This function provides ability to fetch the intron class data for 265
species in MiDB database. For more information refer to the [MiDB
database](https://midb.pnb.uconn.edu/aboutus.php)

## Usage

``` r
searchMidb(org, type = "intron")
```

## Arguments

- org:

  Name of Organism. Can be obtained by loading midbSpecies data from the
  package

- type:

  Type of data to be fetched. Default: "intron"

## Value

df This function returns a dataframe of intron classification in
Plasmodium species.

## Examples

``` r
if (FALSE) { # \dontrun{
load("data/midbSpecies.rda")
## Fetching intron data from MiDB for P. falciparum
df <- searchMidb(midbSpecies$`Available Species`[196])
} # }
```
