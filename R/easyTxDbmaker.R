#' Create a TxDb from VEuPathDB GFF3 annotations
#'
#' Downloads GFF3 and genomic FASTA with an optional VEuPathDB API key,
#' or reads either input from a local file. Temporary downloads are deleted
#' on exit. Returns a TxDb object, not an installable TxDb package.
#'
#' IMPORTANT: `fasta` must contain genomic chromosome/contig sequences, NOT
#' annotated transcript sequences. If `fasta = NULL`, chromosome lengths are
#' not supplied explicitly; the GFF may still contain sequence metadata.
#'
#' @import stringr httr2
#'
#' @param gff Local GFF3/GFF3.gz path or HTTPS URL.
#' @param fasta Optional local genomic FASTA/FASTA.gz path or HTTPS URL.
#'   Do not supply an AnnotatedTranscripts FASTA as chromosome information.
#' @param abbr Genome assembly label stored in the TxDb.
#' @param org Scientific organism name.
#' @param taxid NCBI taxonomy identifier.
#' @param db Description of the annotation source/release.
#' @param api_key Optional VEuPathDB bearer token; defaults to the
#'   VEUPATHDB_API_KEY environment variable. Sent only to trusted VEuPathDB-
#'   family hosts over HTTPS. Public URLs work without a key.
#' @param timeout Positive download timeout in seconds.
#' @param verbose Print status messages.
#'
#' @return A GenomicFeatures TxDb object.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' txdb <- easyTxDbmaker(gff = "https://plasmodb.org/a/service/raw-files/release-71/Pfalciparum3D7/gff/data/PlasmoDB-71_Pfalciparum3D7.gff", fasta = "https://plasmodb.org/a/service/raw-files/release-71/Pfalciparum3D7/fasta/data/PlasmoDB-71_Pfalciparum3D7_Genome.fasta")
#'
#' # GFF-only, without explicit chromosome lengths:
#' txdb <- easyTxDbmaker(gff = "Pfalciparum3D7.gff", fasta = NULL)
#'
#' # Offline use with accurate chromosome lengths:
#' txdb <- easyTxDbmaker(
#'   gff = "Pfalciparum3D7.gff.gz",
#'   fasta = "Pfalciparum3D7_Genome.fasta.gz"
#' )
#' }
easyTxDbmaker <- function(
  gff = NULL,
  fasta = NULL,
  abbr = "Pfalciparum3D7",
  org = "Plasmodium falciparum",
  taxid = 36329L,
  db = "PlasmoDB release 71",
  api_key = Sys.getenv("VEUPATHDB_API_KEY"),
  timeout = 300,
  verbose = TRUE
) {
  pkgs <- c("GenomeInfoDb", "Biostrings", "txdbmaker", "rtracklayer")

  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_pkgs) > 0) {
    stop(
      "Required package(s) missing: ",
      paste(missing_pkgs, collapse = ", "),
      ". Please install them using BiocManager::install().",
      call. = FALSE
    )
  }

  # 1. Validate inputs early.
  is_one_string <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(trimws(x))
  }
  if (!is_one_string(gff)) stop("'gff' must be a nonempty path or URL.", call. = FALSE)
  if (!is.null(fasta) && !is_one_string(fasta)) {
    stop("'fasta' must be NULL or a nonempty genomic FASTA path or URL.", call. = FALSE)
  }
  if (!is_one_string(abbr) || !is_one_string(org) || !is_one_string(db)) {
    stop("'abbr', 'org', and 'db' must be nonempty strings.", call. = FALSE)
  }
  if (length(taxid) != 1L || is.na(taxid) ||
    !grepl("^[0-9]+$", as.character(taxid))) {
    stop("'taxid' must be one NCBI taxonomy ID.", call. = FALSE)
  }
  if (!is_one_string(api_key) && !identical(api_key, "")) {
    stop("'api_key' must be a single string (or empty).", call. = FALSE)
  }
  if (!is.numeric(timeout) || length(timeout) != 1L ||
    is.na(timeout) || !is.finite(timeout) || timeout <= 0) {
    stop("'timeout' must be a positive number of seconds.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
  }

  # The original implementation used AnnotatedTranscripts.fasta to compute
  # chrominfo; transcript lengths do not represent chromosome lengths.
  if (!is.null(fasta) && grepl("AnnotatedTranscripts", fasta, ignore.case = TRUE)) {
    stop("'fasta' points to AnnotatedTranscripts. Supply the matching genomic ",
      "Genome FASTA instead, or set fasta = NULL.",
      call. = FALSE
    )
  }

  # 2. Create a temporary directory and always clean up our own files.
  stage_dir <- tempfile(pattern = "easyTxDb_stage_")
  if (!dir.create(stage_dir)) stop("Cannot create temporary directory.", call. = FALSE)
  on.exit(unlink(stage_dir, recursive = TRUE), add = TRUE)

  # 3. Obtain a URL or local file, preserving a gzip extension for parsers.
  trusted_domains <- c(
    "veupathdb.org", "eupathdb.org", "plasmodb.org", "toxodb.org",
    "tritrypdb.org", "fungidb.org", "giardiadb.org", "amoebadb.org",
    "piroplasmadb.org", "cryptodb.org", "microsporidiadb.org"
  )

  fetch_file <- function(src, prefix, extension) {
    compressed <- grepl("\\.(gz|bgz)([?#].*)?$", src, ignore.case = TRUE)
    dest <- file.path(stage_dir, paste0(prefix, extension, if (compressed) ".gz" else ""))
    is_url <- grepl("^https?://", src, ignore.case = TRUE)

    if (is_url) {
      host <- tolower(sub("^https?://([^/:?#]+).*", "\\1", src, ignore.case = TRUE))
      trusted <- any(vapply(trusted_domains, function(domain) {
        identical(host, domain) || endsWith(host, paste0(".", domain))
      }, logical(1)))
      https <- grepl("^https://", src, ignore.case = TRUE)

      req <- httr2::request(src) |>
        httr2::req_timeout(timeout)

      if (nzchar(api_key) && trusted && https) {
        req <- httr2::req_auth_bearer_token(req, token = api_key)
      } else if (nzchar(api_key) && trusted && !https) {
        warning("API key not sent over unencrypted HTTP.", call. = FALSE)
      } else if (nzchar(api_key) && !trusted && verbose) {
        message("API key not sent to an unrecognised host: ", host)
      }

      if (verbose) message("Downloading ", prefix, ": ", src)
      tryCatch(
        httr2::req_perform(req, path = dest),
        error = function(e) {
          stop("Unable to download ", prefix, " from ", src, ": ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      )
    } else {
      if (!file.exists(src) || dir.exists(src)) {
        stop("Local ", prefix, " file does not exist: ", src, call. = FALSE)
      }
      if (verbose) message("Using local ", prefix, ": ", src)
      if (!file.copy(src, dest, overwrite = TRUE)) {
        stop("Could not copy local ", prefix, ": ", src, call. = FALSE)
      }
    }

    size <- file.info(dest)$size
    if (!file.exists(dest) || is.na(size) || size == 0L) {
      stop("Empty or missing ", prefix, " file: ", src, call. = FALSE)
    }
    if (compressed) {
      con <- file(dest, open = "rb")
      signature <- tryCatch(readBin(con, what = "raw", n = 2L),
        finally = close(con)
      )
      if (!identical(as.integer(signature), c(31L, 139L))) {
        stop("Expected gzip-compressed ", prefix,
          "; received a different format. Check the URL/file.",
          call. = FALSE
        )
      }
    } else {
      first_line <- readLines(dest, n = 1L, warn = FALSE)
      if (length(first_line) &&
        grepl("^\\s*(<!DOCTYPE|<html|\\{\\s*\"(error|message)\")",
          first_line,
          ignore.case = TRUE
        )) {
        stop("The ", prefix, " input appears to be an HTML/JSON error page: ",
          src,
          call. = FALSE
        )
      }
    }
    dest
  }

  gff_file <- fetch_file(gff, "GFF", ".gff3")
  fasta_file <- if (is.null(fasta)) NULL else fetch_file(fasta, "FASTA", ".fasta")

  # 4. Construct chromosome/contig lengths from GENOMIC FASTA only.
  chrominfo <- NULL
  if (!is.null(fasta_file)) {
    if (verbose) message("Reading genomic FASTA and validating sequence names...")
    dna <- Biostrings::readDNAStringSet(fasta_file, format = "fasta")
    fasta_names <- names(dna)
    seq_ids <- sub("[[:space:]].*$", "", fasta_names)

    if (!length(seq_ids) || anyNA(seq_ids) || any(!nzchar(seq_ids)) ||
      anyDuplicated(seq_ids)) {
      stop("FASTA sequence IDs must be nonempty and unique (first header token).",
        call. = FALSE
      )
    }
    seq_lengths <- as.integer(BiocGenerics::width(dna))
    names(seq_lengths) <- seq_ids

    # Ensure the GFF's genomic coordinates use the SAME sequence IDs as FASTA.
    # This detects the mistaken use of a transcript FASTA even when its name
    # does not contain the word 'AnnotatedTranscripts'.
    gr <- rtracklayer::import(gff_file, format = "gff3")
    if (!length(gr)) stop("GFF contains no genomic features.", call. = FALSE)
    gff_chr <- as.character(GenomeInfoDb::seqnames(gr))
    needed <- unique(gff_chr)
    absent <- setdiff(needed, seq_ids)
    if (length(absent)) {
      stop("Genomic FASTA lacks ", length(absent), " GFF sequence IDs (e.g. ",
        paste(head(absent, 6L), collapse = ", "), "). Supply the matching ",
        "genome FASTA, not transcript FASTA.",
        call. = FALSE
      )
    }

    # Detect incorrect genome assemblies or truncated FASTA sequences.
    max_end <- tapply(BiocGenerics::end(gr), gff_chr, max)
    too_short <- names(max_end)[max_end > seq_lengths[names(max_end)]]
    if (length(too_short)) {
      stop("GFF coordinates exceed FASTA lengths for: ",
        paste(head(too_short, 6L), collapse = ", "),
        ". Check that FASTA and GFF use the same assembly/release.",
        call. = FALSE
      )
    }

    chrominfo <- GenomeInfoDb::Seqinfo(
      seqnames = seq_ids,
      seqlengths = unname(seq_lengths),
      genome = abbr
    )
    if (verbose) {
      message(
        "Validated ", length(needed),
        " GFF sequences against the genomic FASTA."
      )
    }
    rm(dna, gr)
  } else if (verbose) {
    message("No genomic FASTA: proceeding without explicit chromosome lengths.")
  }

  # 5. Retain original paths/URLs in metadata, not temporary filenames.
  metadata <- data.frame(
    name = c("GFF/GTF", if (!is.null(fasta)) "Genome FASTA"),
    value = c(gff, if (!is.null(fasta)) fasta),
    stringsAsFactors = FALSE
  )

  if (verbose) message("Building TxDb from GFF3...")
  txdb <- txdbmaker::makeTxDbFromGFF(
    file = gff_file,
    format = "gff3",
    chrominfo = chrominfo,
    metadata = metadata,
    organism = org,
    taxonomyId = as.integer(taxid),
    dataSource = db
  )

  txdb <- GenomeInfoDb::`genome<-`(txdb, value = abbr)
  if (verbose) message("TxDb created successfully.")
  txdb
}
