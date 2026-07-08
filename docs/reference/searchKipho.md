# Fetch Kinases and Phosphatases from Kipho Database

This function provides ability to query your gene IDs to KiPho Database.

## Usage

``` r
searchKipho(org = "pf", type = "kinase")
```

## Arguments

- org:

  Abbreviation of organism of interest.

  - pb: *Plasmodium berghii*

  - pv: *Plasmodium vivax*

  - pf: *Plasmodium falciparum*

  - pc: *Plasmodium chabaudi*

- type:

  Type of protein class i.e. "kinase" or "phosphatase". Default:
  "kinase"

## Value

df This function returns a dataframe of kinases/phosphatases in
Plasmodium species.

## Examples

``` r
if (FALSE) { # \dontrun{
test <- searchKipho(org="pf")
} # }
```
