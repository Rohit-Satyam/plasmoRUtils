# Create a TxDb from VEuPathDB GFF3 annotations

Downloads GFF3 and genomic FASTA with an optional VEuPathDB API key, or
reads either input from a local file. Temporary downloads are deleted on
exit. Returns a TxDb object, not an installable TxDb package.

## Usage

``` r
easyTxDbmaker(
  gff = NULL,
  fasta = NULL,
  abbr = "Pfalciparum3D7",
  org = "Plasmodium falciparum",
  taxid = 36329L,
  db = "PlasmoDB release 71",
  api_key = Sys.getenv("VEUPATHDB_API_KEY"),
  timeout = 300,
  verbose = TRUE
)
```

## Arguments

- gff:

  Local GFF3/GFF3.gz path or HTTPS URL.

- fasta:

  Optional local genomic FASTA/FASTA.gz path or HTTPS URL. Do not supply
  an AnnotatedTranscripts FASTA as chromosome information.

- abbr:

  Genome assembly label stored in the TxDb.

- org:

  Scientific organism name.

- taxid:

  NCBI taxonomy identifier.

- db:

  Description of the annotation source/release.

- api_key:

  Optional VEuPathDB bearer token; defaults to the VEUPATHDB_API_KEY
  environment variable. Sent only to trusted VEuPathDB- family hosts
  over HTTPS. Public URLs work without a key.

- timeout:

  Positive download timeout in seconds.

- verbose:

  Print status messages.

## Value

A GenomicFeatures TxDb object.

## Details

IMPORTANT: `fasta` must contain genomic chromosome/contig sequences, NOT
annotated transcript sequences. If `fasta = NULL`, chromosome lengths
are not supplied explicitly; the GFF may still contain sequence
metadata.

## Examples

``` r
if (FALSE) { # \dontrun{

txdb <- easyTxDbmaker(gff = "https://plasmodb.org/a/service/raw-files/release-71/Pfalciparum3D7/gff/data/PlasmoDB-71_Pfalciparum3D7.gff", fasta = "https://plasmodb.org/a/service/raw-files/release-71/Pfalciparum3D7/fasta/data/PlasmoDB-71_Pfalciparum3D7_Genome.fasta")

# GFF-only, without explicit chromosome lengths:
txdb <- easyTxDbmaker(gff = "Pfalciparum3D7.gff", fasta = NULL)

# Offline use with accurate chromosome lengths:
txdb <- easyTxDbmaker(
  gff = "Pfalciparum3D7.gff.gz",
  fasta = "Pfalciparum3D7_Genome.fasta.gz"
)
} # }
```
