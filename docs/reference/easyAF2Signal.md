# Check Signal peptide fidelity

A convenience function to compute atomic contacts between a predicted
signal peptide and the rest of the protein in AlphaFold2 structures.
This implements the FORTRAN code from Elcock-Lab/AlphaFold2-Signal, with
added features to calculate the median pLDDT and count residues
surpassing the pLDDT threshold for the signal peptide. Unlike the
FORTRAN version, this function highlights cases where low-confidence
signal peptide structures filter all residues, suggesting they may not
represent true outward-facing signal peptides.

## Usage

``` r
easyAF2Signal(pdb, cut_dist = 4, nsignal = 25, bfac_thresh = 90, nskip = 1)
```

## Arguments

- pdb:

  A link to PDB file or path to PDB file if the file is locally present.

- cut_dist:

  Distance cutoff for defining atomic contacts (Default: 4 Angstroms).

- nsignal:

  Number of residues that comprise of N-terminal signal peptide. If
  SignalP prediction is available for your protein, the predicted length
  should be used. (Default: 25).

- bfac_thresh:

  pLDDT threshold for including residues in the atomic contact
  calculation.

- nskip:

  No. of residues immediately next to the cleavage site to exclude from
  atomic contact calculations.

## Value

A data frame, containing statistics about signal peptide and rest of the
protein. Zero atomic or residue-residue contact is indicative of True
positives while non-zero values are putative False Positives.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- easyAF2Signal("https://alphafold.ebi.ac.uk/files/AF-A0A140KXW6-F1-model_v6.pdb")
} # }
```
