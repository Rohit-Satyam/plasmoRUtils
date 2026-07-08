# Get NOISeq::readData ready annotations

A convenience function to make a data frame containing biological
annotations required by NOISeq to run NOISeq::readData function using
custom GTF/GFF file and FASTA.

## Usage

``` r
easyNOISeqAnnot(
  gff,
  fasta,
  name = "Tgondii",
  select = c("protein_coding_gene", "ncRNA_gene", "pseudogene"),
  geneidcol = "ID",
  genetype = "ebi_biotype"
)
```

## Arguments

- gff:

  A link to GFF file from VEuPathDB or path to GTF file produced using
  AGAT.

- fasta:

  A link or path to the genome fasta file.

- name:

  Name of the organism.

- select:

  Type of features to be selected. By default we select
  "protein_coding_gene","ncRNA_gene"and "pseudogene" as they cover all
  the genes in the VEuPathDB annotation files. If using AGAT formatted
  GTF file, using "gene" would be sufficient as AGAT put all the genes
  types under gene tag.

- geneidcol:

  Use the tag that refers to the gene IDs in the GTF/GFF file.

- genetype:

  Use the tag that refers to the gene subtypes such as "ebi_biotype".

## Value

A dataframe containing annotations per gene such as GC content, gene
description, gene-length, gene start, gene end coordinates and
chromosome information.

## Examples

``` r
if (FALSE) { # \dontrun{
 df <- easyNOISeqAnnot(
 gff="https://toxodb.org/common/downloads/release-68/EpraecoxHoughton/gff/data/ToxoDB-68_EpraecoxHoughton.gff",
 fasta = "https://toxodb.org/common/downloads/release-68/EpraecoxHoughton/fasta/data/ToxoDB-68_EpraecoxHoughton_Genome.fasta")


} # }
```
