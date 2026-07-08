# Quickly make TxDb object

A convenience function to make OrgDb packages quickly to be used with
other GO enrichment packages such as ClusterProfiler.

## Usage

``` r
easyTxDbmaker(
  gff,
  fasta,
  abbr = "TgondiiME49",
  org = "Toxoplasma gondii",
  taxid = 508771,
  db = "ToxoDB release 68"
)
```

## Arguments

- gff:

  A link to GFF file from VEuPathDB or path to GFF file if the file is
  present locally.

- fasta:

  Genome FASTA file obtained from VEuPathDB.

- abbr:

  Abbreviation of the organism.

- org:

  Name of the organism including genus and species.

- taxid:

  NCBI Taxonomy ID. This can be obtained from
  https://www.ncbi.nlm.nih.gov/taxonomy.

- db:

  Name of the database with release information.

## Value

A tar.gz file that can be installed as a package and can be used with GO
enrichment tools such as ClusterProfiler.

## Examples

``` r
if (FALSE) { # \dontrun{
 txdb<-easyTxDbmaker(
 gff="https://toxodb.org/common/downloads/release-68/TgondiiME49/gff/data/ToxoDB-68_TgondiiME49.gff",
 fasta="https://toxodb.org/common/downloads/release-68/TgondiiME49/fasta/data/ToxoDB-68_TgondiiME49_Genome.fasta",
 abbr="TgondiiME49",
 taxid=508771,org = "Toxoplasma gondii ME49",
 db = "ToxoDB release 68")


} # }
```
