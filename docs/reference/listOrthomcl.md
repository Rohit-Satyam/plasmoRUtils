# List species metadata present in OrthoMCL

A convenience function to quickly fetch the species related vocabulary
of 713 species used by [OrthoMCL](https://orthomcl.org/orthomcl/app)
database. This function helps users choose IDs of the organisms between
which they wish to fetch the paired orthologs. See also:
[`getpairedOrthologs()`](https://rohit-satyam.github.io/plasmoRUtils/reference/getpairedOrthologs.md).

## Usage

``` r
listOrthomcl()
```

## Value

A data frame containing information about species present in
InParanoiDB9 and their taxon ID.

## Examples

``` r
if (FALSE) { # \dontrun{
listOrthomcl()
} # }
```
