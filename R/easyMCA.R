#' Download data from the Malaria Cell Atlas
#'
#' Download an MCA dataset archive, or read a locally available ZIP archive,
#' and import the requested metadata or expression table.
#'
#' @param url URL of an MCA ZIP archive (for example, one returned by
#'   [listMCA()]), or the path to a local ZIP archive.
#' @param type Character string specifying the table to import: `"data"`
#'   (cell metadata; default), `"exp"` (processed expression values), or
#'   `"raw"` (raw counts).
#' @param timeout Download timeout in seconds (default: 300). This applies
#'   only to remote URLs; the previous R timeout option is restored on exit.
#'
#' @return A `data.frame` containing the selected table, with the first CSV
#'   column used as row names. For expression tables, genes are expected to
#'   occupy rows and cells columns, as provided by the MCA download.
#'
#' @details Only the requested CSV file is extracted, into an isolated
#'   temporary directory. Temporary files are removed when the function exits,
#'   including if an error occurs. Local input ZIP archives are never deleted.
#'   The supported archive members end in `-data.csv`, `-exp.csv`, or
#'   `-raw.csv` (also accepting underscore separators and `.csv.gz`).
#'
#' @md
#' @export
#'
#' @examples
#' \dontrun{
#' url <- "https://www.malariacellatlas.org/downloads/pf-ch10x-set4-biorxiv.zip"
#'
#' metadata <- easyMCA(url, type = "data")
#' expression <- easyMCA(url, type = "exp")
#' raw_counts <- easyMCA(url, type = "raw")
#'
#' # An already downloaded archive can also be used offline.
#' metadata <- easyMCA("path/to/pf-ch10x-set4-biorxiv.zip", "data")
#'
#' # When creating a Seurat object, check that cells match and
#' # convert the counts data frame to a numeric matrix.
#' stopifnot(all(colnames(raw_counts) %in% rownames(metadata)))
#' seurat_obj <- Seurat::CreateSeuratObject(
#'   counts = as.matrix(raw_counts),
#'   meta.data = metadata[colnames(raw_counts), , drop = FALSE],
#'   project = "MCA"
#' )
#' }
easyMCA <- function(url, type = "data", timeout = 300) {
  # 1. Validate arguments early, before creating or downloading any files.
  if (!is.character(url) || length(url) != 1L ||
    is.na(url) || !nzchar(trimws(url))) {
    stop("'url' must be a single non-empty URL or local ZIP path.",
      call. = FALSE
    )
  }

  choices <- c("data", "exp", "raw")
  if (!is.character(type) || length(type) != 1L ||
    is.na(type) || !type %in% choices) {
    stop("'type' must be exactly one of: ",
      paste(shQuote(choices), collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (!is.numeric(timeout) || length(timeout) != 1L ||
    is.na(timeout) || !is.finite(timeout) || timeout <= 0) {
    stop("'timeout' must be a single positive number of seconds.",
      call. = FALSE
    )
  }

  # 2. Work only in an isolated directory. In particular, do not delete
  #    the original ZIP if the user supplies a local path.
  work_dir <- tempfile(pattern = "easyMCA_")
  if (!dir.create(work_dir)) {
    stop("Could not create a temporary directory.", call. = FALSE)
  }
  on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

  is_url <- grepl("^https?://", url, ignore.case = TRUE)

  if (is_url) {
    zip_file <- file.path(work_dir, "mca_download.zip")

    old_timeout <- getOption("timeout")
    options(timeout = max(old_timeout, timeout))
    on.exit(options(timeout = old_timeout), add = TRUE)

    status <- tryCatch(
      utils::download.file(url,
        destfile = zip_file,
        mode = "wb", quiet = TRUE
      ),
      error = function(e) {
        stop("Unable to download MCA archive: ", conditionMessage(e),
          call. = FALSE
        )
      }
    )

    if (!identical(as.integer(status), 0L)) {
      stop("MCA download failed (status: ", status, ").", call. = FALSE)
    }
  } else {
    if (grepl("^[A-Za-z][A-Za-z0-9+.-]*://", url)) {
      stop("Only HTTP(S) URLs and local ZIP paths are supported.",
        call. = FALSE
      )
    }
    if (!file.exists(url) || dir.exists(url)) {
      stop("Local MCA archive not found: ", url, call. = FALSE)
    }
    zip_file <- normalizePath(url, mustWork = TRUE)
  }

  size <- file.info(zip_file)$size
  if (is.na(size) || size == 0) {
    stop("The MCA ZIP archive is empty or unavailable.", call. = FALSE)
  }

  # 3. Inspect ZIP contents WITHOUT extracting everything.
  members <- tryCatch(
    utils::unzip(zip_file, list = TRUE),
    error = function(e) {
      stop("Cannot read MCA ZIP archive: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (!is.data.frame(members) || !"Name" %in% names(members) ||
    nrow(members) == 0L) {
    stop("The ZIP archive contains no readable members.", call. = FALSE)
  }

  member_names <- as.character(members$Name)

  # Match the requested type as the FINAL filename token; e.g.,
  # dataset-raw.csv matches 'raw', but dataset-raw-data.csv does not.
  pattern <- paste0("(^|[-_])", type, "\\.csv(\\.gz)?$")
  candidates <- grepl(pattern, basename(member_names),
    ignore.case = TRUE, perl = TRUE
  )

  # Exclude operating-system metadata files.
  candidates <- candidates &
    !grepl("(^|/)__MACOSX(/|$)|(^|/)\\._", member_names)

  matches <- member_names[candidates]

  if (length(matches) == 0L) {
    stop("No CSV file for type '", type, "' found in the MCA archive. ",
      "Available entries: ",
      paste(utils::head(member_names, 20L), collapse = ", "),
      call. = FALSE
    )
  }
  if (length(matches) > 1L) {
    stop("Multiple CSV files match type '", type, "': ",
      paste(matches, collapse = ", "),
      ". Please use an unambiguous archive.",
      call. = FALSE
    )
  }

  target <- matches[[1L]]

  # 4. Reject unsafe ZIP paths (absolute paths, drive paths and '..').
  if (grepl("\\\\", target) ||
    grepl("^(/|[A-Za-z]:)", target) ||
    grepl("(^|/)\\.\\.(/|$)", target) ||
    grepl("(^|/)\\.(/|$)", target)) {
    stop("Unsafe member path in ZIP archive: ", target, call. = FALSE)
  }

  extract_dir <- file.path(work_dir, "extracted")
  if (!dir.create(extract_dir)) {
    stop("Could not create extraction directory.", call. = FALSE)
  }

  # Extract only the selected file, not the whole dataset.
  tryCatch(
    utils::unzip(zip_file,
      files = target, exdir = extract_dir,
      overwrite = FALSE, unzip = "internal"
    ),
    error = function(e) {
      stop("Could not extract '", target, "': ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  selected_file <- file.path(extract_dir, target)
  if (!file.exists(selected_file) || dir.exists(selected_file)) {
    stop("Selected CSV was not extracted: ", target, call. = FALSE)
  }

  # Check the actual resolved path too, in case an archive contains links.
  root <- normalizePath(extract_dir, winslash = "/", mustWork = TRUE)
  selected_file <- normalizePath(selected_file,
    winslash = "/",
    mustWork = TRUE
  )
  if (!startsWith(selected_file, paste0(root, "/"))) {
    stop("Extracted member resolves outside the temporary directory.",
      call. = FALSE
    )
  }

  # 5. Read the selected CSV. A gzipped member is read through a connection,
  #    so no additional decompression package is required.
  if (grepl("\\.csv\\.gz$", target, ignore.case = TRUE)) {
    con <- gzfile(selected_file, open = "rt")
    result <- tryCatch(
      utils::read.csv(con,
        header = TRUE, check.names = FALSE,
        stringsAsFactors = FALSE
      ),
      error = function(e) {
        stop("Failed to read MCA CSV: ", conditionMessage(e),
          call. = FALSE
        )
      },
      finally = close(con)
    )
  } else {
    result <- tryCatch(
      data.table::fread(selected_file,
        sep = ",", header = TRUE,
        data.table = FALSE, check.names = FALSE,
        showProgress = FALSE
      ),
      error = function(e) {
        stop("Failed to read MCA CSV: ", conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }

  # 6. Convert the first column to row names without changing the remaining
  #    metadata or numeric expression columns.
  if (nrow(result) == 0L || ncol(result) < 2L) {
    stop("The selected MCA table has no rows or too few columns.",
      call. = FALSE
    )
  }

  row_ids <- as.character(result[[1L]])
  if (anyNA(row_ids) || any(!nzchar(trimws(row_ids)))) {
    stop("The first CSV column contains missing or blank row identifiers.",
      call. = FALSE
    )
  }
  if (anyDuplicated(row_ids)) {
    stop("The first CSV column contains duplicate row identifiers; ",
      "cannot safely assign row names.",
      call. = FALSE
    )
  }

  result <- result[-1L]
  rownames(result) <- row_ids

  if (type %in% c("exp", "raw") &&
    any(!vapply(result, is.numeric, logical(1)))) {
    warning("Some expression columns are not numeric; inspect the MCA CSV ",
      "before downstream analysis.",
      call. = FALSE
    )
  }

  return(result)
}
