# Quick mapping of PDB IDs to Uiprot IDs

A convenience function to quickly convert the PDB IDs to Uniprot IDs. If
a protein is multimeric, corresponding Uniprot IDs are returned for
them. This function uses PDBe API.

## Usage

``` r
pdb2uniprot(pdbid, timeout = 60, max_tries = 3)
```

## Arguments

- pdbid:

  A single PDB Id.

- timeout:

  Timeout time in seconds

- max_tries:

  Maximum retries allowed.

## Value

A data frame, Uniprot IDs, Chain IDs and start and end coordinates of
the chains.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- pdb2uniprot("9FIA")
} # }
```
