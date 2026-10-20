#' Quickly create an organism-specific OrgDb package
#'
#' Create an AnnotationDbi-compatible OrgDb source package from a VEuPathDB
#' GFF3 file and a Gene Ontology Association File (GAF). Both inputs may be
#' local files or HTTP(S) URLs; gzip-compressed inputs are supported.
#'
#' The GFF-derived genes define the annotation universe. Negative (`NOT`) GO
#' annotations, invalid/unknown GO identifiers, and associations whose gene
#' IDs do not match the GFF are excluded. The full GFF gene universe is kept,
#' including genes without a GO annotation.
#'
#' Note: If the function reports too many GO IDs as obsolete, try using the latest GAF file from VEuPathDB. Obsolete terms are routinely removed i.e. version 68 will have more obsolete GO IDs and associations than version 71. Also, routinely update GO.db and Annotationdbi R packages
#'
#'
#' @param gff Local file path or HTTP(S) URL to a GFF3 file (.gff/.gff3,
#'   optionally gzip-compressed).
#' @param gaf Local file path or HTTP(S) URL to a GAF file (plain or .gz).
#' @param out.dir Directory where the package source and tar.gz will be saved.
#' @param taxid NCBI taxonomy identifier.
#' @param genus Organism genus, used in the generated package name.
#' @param sp Species (optionally including strain), without whitespace.
#' @param version Valid R package version, e.g., "0.1.0".
#' @param verbose Whether to print informational messages.
#' @param maintainer Name and email, e.g., "Jane Doe <jane@example.org>".
#' @param api_key Optional VEuPathDB bearer token. It is sent only to
#'   recognised VEuPathDB-family HTTPS hosts. The requested endpoint must
#'   support bearer-token authentication.
#' @param timeout Positive download timeout, in seconds.
#' @param overwrite Whether to replace an existing package of the same name
#'   and version. Defaults to FALSE.
#'
#' @return Invisibly returns the normalized path of the generated source
#'   package archive (.tar.gz). The unpacked package source is also retained.
#'
#' @examples
#' \dontrun{
#' archive <- easymakeOrgdb(
#'   gff = "PlasmoDB-68_Pfalciparum3D7.gff",
#'   gaf = "PlasmoDB-68_Pfalciparum3D7_Curated_GO.gaf.gz",
#'   out.dir = "./annotation",
#'   taxid = 36329,
#'   genus = "Plasmodium",
#'   sp = "falciparum3D7",
#'   version = "0.1.0",
#'   maintainer = "Jane Doe <jane@example.org>",
#'   verbose = TRUE, overwrite = TRUE
#' )
#' install.packages(archive, repos = NULL, type = "source")
#' }
#'
#' @export
easymakeOrgdb <- function(
  gff,
  gaf,
  out.dir = ".",
  taxid = 36329,
  genus = "Plasmodium",
  sp = "falciparum3D7",
  version = "0.1.0",
  verbose = TRUE,
  maintainer = "John Doe <johndoe@gmail.com>",
  api_key = Sys.getenv("VEUPATHDB_API_KEY"),
  timeout = 300,
  overwrite = FALSE
) {
  # 1. Validate arguments --------------------------------------------------
  scalar_string <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }

  if (!scalar_string(gff) || !scalar_string(gaf)) {
    stop("'gff' and 'gaf' must each be a nonempty file path or URL.",
      call. = FALSE
    )
  }
  if (!scalar_string(out.dir) || !scalar_string(genus) ||
    !scalar_string(sp) || !scalar_string(maintainer)) {
    stop("'out.dir', 'genus', 'sp', and 'maintainer' must be nonempty strings.",
      call. = FALSE
    )
  }
  if (!is.numeric(timeout) || length(timeout) != 1L ||
    is.na(timeout) || !is.finite(timeout) || timeout <= 0) {
    stop("'timeout' must be one positive, finite number.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose) ||
    !is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("'verbose' and 'overwrite' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!scalar_string(api_key) &&
    !(is.character(api_key) && length(api_key) == 1L &&
      !is.na(api_key) && identical(api_key, ""))) {
    stop("'api_key' must be a single string (empty is allowed).", call. = FALSE)
  }

  version <- as.character(version)
  if (!scalar_string(version)) {
    stop("'version' must be a single valid package-version string.",
      call. = FALSE
    )
  }
  tryCatch(
    base::package_version(version),
    error = function(e) {
      stop("Invalid package version '", version, "': ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  # Species is used in AnnotationForge's generated package name.
  if (!grepl("^[A-Za-z][A-Za-z0-9]*$", genus) ||
    !grepl("^[A-Za-z][A-Za-z0-9.]*$", sp)) {
    stop("'genus' and 'sp' must contain valid package-name components ",
      "(without spaces, hyphens or underscores).",
      call. = FALSE
    )
  }
  taxid <- as.character(taxid)
  if (!scalar_string(taxid) || !grepl("^[0-9]+$", taxid) ||
    is.na(suppressWarnings(as.integer(taxid))) || as.integer(taxid) <= 0L) {
    stop("'taxid' must be a positive NCBI taxonomy identifier.",
      call. = FALSE
    )
  }

  if (!dir.exists(out.dir)) {
    if (!dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create output directory: ", out.dir, call. = FALSE)
    }
  }
  out.dir <- normalizePath(out.dir, mustWork = TRUE)
  if (file.access(out.dir, mode = 2L) != 0L) {
    stop("Output directory is not writable: ", out.dir, call. = FALSE)
  }

  # 2. Fetch files ---------------------------------------------------------
  # A temporary staging area also keeps incomplete builds out of out.dir.
  stage_dir <- tempfile(pattern = ".orgdb_stage_", tmpdir = out.dir)
  if (!dir.create(stage_dir)) {
    stop("Cannot create temporary build directory in: ", out.dir)
  }
  on.exit(unlink(stage_dir, recursive = TRUE), add = TRUE)

  trusted_domains <- c(
    "veupathdb.org", "eupathdb.org", "plasmodb.org", "toxodb.org",
    "tritrypdb.org", "fungidb.org", "giardiadb.org", "amoebadb.org",
    "piroplasmadb.org", "cryptodb.org", "microsporidiadb.org"
  )

  fetch_file <- function(src, dest) {
    remote <- grepl("^https?://", src, ignore.case = TRUE)
    if (remote) {
      req <- httr2::request(src) |>
        httr2::req_timeout(timeout)
      host <- tolower(sub("^https?://([^/:?#]+).*", "\\1", src,
        ignore.case = TRUE
      ))
      trusted <- any(vapply(trusted_domains, function(domain) {
        identical(host, domain) || endsWith(host, paste0(".", domain))
      }, logical(1)))
      if (nzchar(api_key) &&
        grepl("^https://", src, ignore.case = TRUE) && trusted) {
        req <- httr2::req_auth_bearer_token(req, token = api_key)
      }
      if (verbose) message("Downloading: ", src)
      httr2::req_perform(req, path = dest)
    } else {
      if (!file.exists(src) || dir.exists(src)) {
        stop("Input file does not exist: ", src, call. = FALSE)
      }
      if (!file.copy(src, dest, overwrite = TRUE)) {
        stop("Could not copy local input: ", src, call. = FALSE)
      }
    }
    info <- file.info(dest)
    if (is.na(info$size) || info$size == 0L) {
      stop("Downloaded/copied input is empty: ", src, call. = FALSE)
    }
    dest
  }

  # Detect gzip using its file signature rather than trusting URL extensions.
  # This covers URLs with query strings and locally misnamed compressed files.
  is_gzip <- function(path) {
    con <- file(path, open = "rb")
    on.exit(close(con), add = TRUE)
    identical(as.integer(readBin(con, what = "raw", n = 2L)), c(31L, 139L))
  }

  gff_raw <- fetch_file(gff, file.path(stage_dir, "gff_download"))
  gaf_raw <- fetch_file(gaf, file.path(stage_dir, "gaf_download"))
  gff_tmp <- file.path(stage_dir, if (is_gzip(gff_raw)) "input.gff3.gz" else "input.gff3")
  gaf_tmp <- file.path(stage_dir, if (is_gzip(gaf_raw)) "input.gaf.gz" else "input.gaf")
  if (!file.rename(gff_raw, gff_tmp) || !file.rename(gaf_raw, gaf_tmp)) {
    stop("Could not prepare temporary annotation files.", call. = FALSE)
  }

  # 3. Parse GFF3 ----------------------------------------------------------
  if (verbose) message("Importing GFF3 annotations...")
  gff_data <- rtracklayer::import.gff3(gff_tmp)
  feature_type <- as.character(gff_data$type)

  if ("protein_coding_gene" %in% feature_type) {
    gene_features <- gff_data[feature_type == "protein_coding_gene"]
  } else if ("gene" %in% feature_type) {
    gene_features <- gff_data[feature_type == "gene"]
    attributes <- names(S4Vectors::mcols(gene_features))
    biotype_col <- intersect(c("gene_biotype", "biotype", "gene_type"), attributes)
    if (length(biotype_col)) {
      biotype <- as.character(S4Vectors::mcols(gene_features)[[biotype_col[1L]]])
      gene_features <- gene_features[
        !is.na(biotype) & biotype %in% c("protein_coding", "protein_coding_gene")
      ]
    } else {
      warning("GFF gene records lack a recognised biotype attribute; ",
        "including all 'gene' features.",
        call. = FALSE
      )
    }
  } else {
    stop("No 'protein_coding_gene' or 'gene' features found in GFF.",
      call. = FALSE
    )
  }
  if (!length(gene_features)) {
    stop("No genes remain after GFF feature filtering.", call. = FALSE)
  }
  attributes <- names(S4Vectors::mcols(gene_features))
  if (!"ID" %in% attributes) {
    stop("GFF gene features lack an ID attribute.", call. = FALSE)
  }
  symbols <- if ("Name" %in% attributes) {
    as.character(gene_features$Name)
  } else {
    rep(NA_character_, length(gene_features))
  }

  gene_info <- data.frame(
    GID = as.character(gene_features$ID),
    SYMBOL = symbols,
    CHR = as.character(GenomeInfoDb::seqnames(gene_features)),
    START = BiocGenerics::start(gene_features),
    END = BiocGenerics::end(gene_features),
    STRAND = as.character(GenomicRanges::strand(gene_features)),
    stringsAsFactors = FALSE
  )
  valid_ids <- !is.na(gene_info$GID) & nzchar(trimws(gene_info$GID))
  gene_info <- gene_info[valid_ids, , drop = FALSE]
  if (!nrow(gene_info)) stop("No valid gene IDs in GFF.", call. = FALSE)
  missing_symbols <- is.na(gene_info$SYMBOL) | !nzchar(trimws(gene_info$SYMBOL))
  gene_info$SYMBOL[missing_symbols] <- gene_info$GID[missing_symbols]

  # AnnotationForge requires unique table rows; a conflicting duplicate gene
  # ID is ambiguous and must be resolved in the input annotation file.
  gene_info <- unique(gene_info)
  if (anyDuplicated(gene_info$GID)) {
    bad <- unique(gene_info$GID[duplicated(gene_info$GID)])
    stop("GFF has conflicting records for gene ID(s): ",
      paste(head(bad, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  if (anyNA(gene_info[, c("CHR", "START", "END", "STRAND")])) {
    stop("Some GFF genes have missing chromosome/coordinate/strand data.",
      call. = FALSE
    )
  }

  # 4. Parse and curate GAF ------------------------------------------------
  if (verbose) message("Importing GAF annotations...")
  gaf_data <- readr::read_tsv(
    gaf_tmp,
    comment = "!", col_names = FALSE,
    col_types = readr::cols(.default = readr::col_character()),
    trim_ws = TRUE, progress = FALSE, show_col_types = FALSE
  )
  if (ncol(gaf_data) < 15L || nrow(gaf_data) == 0L) {
    stop("Invalid/empty GAF: at least 15 tab-separated columns are required.",
      call. = FALSE
    )
  }
  # Standard GAF: column 2 = DB Object ID; 4 = qualifier/relation;
  # 5 = GO ID; 7 = evidence code. If a gene product is highly similar to a protein family but has been experimentally proven not to have the expected function or location, curators use the NOT qualifier.
  qualifier <- trimws(gaf_data[[4L]])
  qualifier[is.na(qualifier)] <- ""
  negative <- grepl("(^|\\|)NOT(\\||$)", qualifier, perl = TRUE)
  n_negative <- sum(negative)
  gaf_data <- gaf_data[!negative, , drop = FALSE]

  go_mappings <- data.frame(
    GID = trimws(gaf_data[[2L]]),
    GO = trimws(gaf_data[[5L]]),
    EVIDENCE = trimws(gaf_data[[7L]]),
    stringsAsFactors = FALSE
  )
  complete_go <- !is.na(go_mappings$GID) & nzchar(go_mappings$GID) &
    !is.na(go_mappings$GO) & grepl("^GO:[0-9]{7}$", go_mappings$GO) &
    !is.na(go_mappings$EVIDENCE) & nzchar(go_mappings$EVIDENCE)
  go_mappings <- go_mappings[complete_go, , drop = FALSE]
  if (!nrow(go_mappings)) {
    stop("No valid GO associations remain after GAF filtering.", call. = FALSE)
  }

  gaf_ids <- unique(go_mappings$GID)
  matching_ids <- intersect(gaf_ids, gene_info$GID)
  if (!length(matching_ids)) {
    stop("No GAF IDs match the GFF gene IDs. Check releases and ID systems.",
      call. = FALSE
    )
  }
  if (length(matching_ids) / length(gaf_ids) < 0.5) {
    warning("Fewer than 50% of GAF gene IDs match the GFF. ",
      "Check species, strains, releases and ID types.",
      call. = FALSE
    )
  }
  # Keep every GFF gene in gene_info but only matching GO associations.
  go_mappings <- go_mappings[go_mappings$GID %in% matching_ids, , drop = FALSE]

  # AnnotationForge expands GO terms using installed GO.db; unknown GOIDs
  # can cause ontology-mapping failures before its own GO-term filtering step.
  recognised_ids <- AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID")
  recognised <- go_mappings$GO %in% recognised_ids
  missing <- go_mappings[!go_mappings$GO %in% recognised_ids, ]
  n_unrecognised <- nrow(missing)
  if (n_unrecognised) {
    warning(n_unrecognised, " associations use GO IDs absent from the ",
      "installed GO.db and were removed.",
      call. = FALSE
    )
  }
  go_mappings <- unique(go_mappings[recognised, , drop = FALSE])
  if (!nrow(go_mappings)) {
    stop("No recognised GO associations remain.", call. = FALSE)
  }

  if (verbose) {
    message("GFF genes: ", nrow(gene_info))
    message("GFF genes with GO: ", length(unique(go_mappings$GID)))
    message("GO associations: ", nrow(go_mappings))
    message("Negative associations excluded: ", n_negative)
    message("Unrecognised GO associations excluded: ", n_unrecognised)
    message("Unique unrecognised GO IDs:", length(unique(missing$GO)))
    message(
      "Unrecognised GO associations by evidence code:\n",
      paste(
        sprintf(
          "  %-5s : %d",
          names(table(missing$EVIDENCE)),
          as.integer(table(missing$EVIDENCE))
        ),
        collapse = "\n"
      )
    )
  }

  # 5. Build in staging directory -----------------------------------------
  # makeOrgPackage() returns the *actual* generated source directory path.
  # This avoids fragile assumptions about AnnotationForge's name generator.
  pkg_dir <- AnnotationForge::makeOrgPackage(
    chromosome = gene_info,
    go = go_mappings,
    version = version,
    maintainer = maintainer,
    author = maintainer,
    outputDir = stage_dir,
    tax_id = taxid,
    genus = genus,
    species = sp,
    goTable = "go",
    verbose = verbose
  )
  if (!dir.exists(pkg_dir)) {
    stop("AnnotationForge failed to create the package source directory.",
      call. = FALSE
    )
  }

  # pkgbuild invokes R CMD build (no Unix-specific tar/rm shell commands).
  staged_archive <- pkgbuild::build(
    path = pkg_dir,
    dest_path = stage_dir,
    quiet = !verbose
  )
  if (!file.exists(staged_archive)) {
    stop("R source-package build did not produce a tar.gz file.", call. = FALSE)
  }

  # Only touch existing destination files AFTER a successful staged build.
  final_dir <- file.path(out.dir, basename(pkg_dir))
  final_archive <- file.path(out.dir, basename(staged_archive))
  existing <- file.exists(final_dir) || file.exists(final_archive)
  if (existing && !overwrite) {
    stop("Output package already exists: ", final_dir,
      ". Use overwrite = TRUE to replace it.",
      call. = FALSE
    )
  }
  if (overwrite) {
    # Save previous outputs in staging for rollback if moving new outputs fails.
    backup_dir <- file.path(stage_dir, "previous_source")
    backup_archive <- file.path(stage_dir, "previous_archive.tar.gz")
    had_dir <- file.exists(final_dir)
    had_archive <- file.exists(final_archive)
    if (had_dir && !file.rename(final_dir, backup_dir)) {
      stop("Could not back up existing source directory.", call. = FALSE)
    }
    if (had_archive && !file.rename(final_archive, backup_archive)) {
      if (had_dir) file.rename(backup_dir, final_dir)
      stop("Could not back up existing archive.", call. = FALSE)
    }
  } else {
    had_dir <- had_archive <- FALSE
  }

  copied_dir <- file.rename(pkg_dir, final_dir)
  copied_archive <- if (copied_dir) file.rename(staged_archive, final_archive) else FALSE
  if (!copied_dir || !copied_archive) {
    # Roll back the transaction as far as possible.
    if (copied_dir && file.exists(final_dir)) {
      file.rename(final_dir, pkg_dir)
    }
    if (copied_archive && file.exists(final_archive)) {
      file.rename(final_archive, staged_archive)
    }
    restored <- TRUE
    if (had_dir) restored <- file.rename(backup_dir, final_dir) && restored
    if (had_archive) {
      restored <- file.rename(backup_archive, final_archive) && restored
    }
    if (!restored) {
      # Preserve backups if restoration failed; don't delete the staging dir.
      on.exit(NULL, add = FALSE)
      warning("Rollback was incomplete. Backups are in: ", stage_dir,
        call. = FALSE
      )
    }
    stop("Could not place completed OrgDb package in output directory.",
      call. = FALSE
    )
  }

  final_archive <- normalizePath(final_archive, mustWork = TRUE)
  if (verbose) message("OrgDb package created: ", final_archive)
  invisible(final_archive)
}
