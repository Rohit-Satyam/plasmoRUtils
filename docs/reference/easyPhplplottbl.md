# Quickly plot PhenoPlasm summary tables

This function generates Disruptability and Mutant Phenotype tables in R,
mirroring the style of Phenoplasm visualizations.

## Usage

``` r
easyPhplplottbl(file, skip = 2)
```

## Arguments

- file:

  Path of Phenotype.txt file obtained from PhenoPlasm database or a data
  frame.

- skip:

  Number of lines to skip in the file. Default 2.

## Value

A gt table plot.

## Examples

``` r
if (FALSE) { # \dontrun{

## Read the table generated from Phenoplasm advance search and pass the resulting data frame
df <- read.csv("phenotype.txt", skip = 2, sep = "\t") %>%
dplyr::select(-3, -4) %>%
dplyr::rename_with(~ gsub("Sprozoite", "Sporozoite", .x))
easyPhplplottbl(df)

## Pass the file path directly
easyPhplplottbl("phenotype.txt")

## Load example data frame
data(pf3d7PhplTable)

easyPhplplottbl(pf3d7PhplTable)
} # }
```
