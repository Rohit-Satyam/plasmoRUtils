#' Performing quick ORA analysis using Fisher or KS and an explicit background
#'
#' Runs topGO enrichment using an explicit eligible-gene universe. Fisher tests
#' over-representation of selected genes; KS uses genuine gene-level p-value
#' ranks for the complete background. Results can be plotted using
#' [easyGOPlot()].
#'
#' @import dplyr tibble
#' @importFrom readr read_tsv
#'
#' @param geneID Uniquely named numeric vector of p-values in, named by
#'   gene ID. For Fisher, either selected genes or the complete tested vector
#'   may be supplied; genes with p-value < `cutoff` are selected. For KS, supply
#'   the actual p-value for EVERY gene in `bkggset` (not just DEGs).
#' @param bkggset Nonempty character vector of all genes eligible for the
#'   analysis, typically genes eligible under your RNA-seq selection procedure.
#'   Required for both tests.
#' @param gaf Path or URL to a GO Annotation File (GAF; .gz is supported).
#'   GAF column 2 must contain gene IDs corresponding to `bkggset`. When using this argument, set useBiomart=FALSE and useGAF=TRUE.
#' @param useBiomart Logical To enable usage of BiomaRt to fetch GO terms. Default: TRUE.
#' @param useGAF Logical To enable usage of custom .gaf file to fetch GO terms. Default: FALSE
#' @param OrgDb An OrgDb object loaded by the caller, e.g. after running
#'   `library(org.Pfalciparum3D7.eg.db)`. Passing it selects offline GO
#'   annotations instead of the default BioMart route.
#' @param orgdb_keytype The key type of the IDs in `bkggset` and `geneID`.
#'   Default `"GID"` for custom packages built with AnnotationForge; check
#'   the keytypes available in your particular database.
#' @param mart,gset BioMart name and dataset, respectively.
#' @param algo Argument to specify which algorithm to be used for enrichment by topGO. For possible options use [topGO::whichAlgorithms()]
#' @param stats Argument to specify the statistical test to be used for enrichment by topGO. For possible values use [topGO::whichTests()]. Default: "fisher".
#' @param category Specify the category of Over-representation analysis such as "BP" for Biological Process, "MF" for Molecular Function and "CC" for Cellular Component Enrichment.
#' @param fdr logical. Perform multiple testing correction testing. Default (FALSE)
#' @param correction Adjustment method accepted by [stats::p.adjust()]. Possible values: ""holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr". Read section 6.2 in topGO documentation before performing correction. Correction when using elim and weight is usually not recommended.
#' @param cutoff P-value threshold defining selected genes, default 0.05.
#'   For KS this changes selected-gene reporting, not the KS statistic.
#' @param go_cutoff Threshold for reporting raw GO p-values or adjusted
#'   GO p-values when `fdr = TRUE`.
#' @param return_GOdata If TRUE, return the underlying `topGOdata` object.
#'
#' @details The OrgDb route uses GO (direct annotations), leaving ancestor
#'   propagation to topGO. topGO itself depends on AnnotationDbi, so no
#'   additional installation or declared direct Import is needed for this route.
#' @return A data.frame with GO IDs, terms, numeric `fisher` or `ks` results,
#'   optional `padj`, and associated/selected genes. The `gene_universe`
#'   attribute gives available, feasible and selected-feasible counts. When
#'   `return_GOdata = TRUE`, returns a `topGOdata` object. Full descriptions
#'   of GO terms truncated by GenTable are restored when AnnotationDbi is installed.
#' @export
#' @examples
#' \dontrun{
#' # res: DESeq2 results with Geneid, pvalue, padj and log2FoldChange.
#' ok <- !is.na(res$padj) & !is.na(res$Geneid)
#' background.gset <- res$Geneid[ok]
#' down <- ok & !is.na(res$log2FoldChange) &
#'   res$padj < 0.05 & res$log2FoldChange < 0
#' geneList <- stats::setNames(res$padj[down], res$Geneid[down])
#'
#' # Fisher accepts a DEG-only list; the full universe is background.gset.
#' gores.fisher <- easytopGO(
#'   geneID = geneList, bkggset = background.gset,
#'   gaf = "annotation.gaf.gz", useBiomart = FALSE, useGAF = TRUE,
#'   stats = "fisher", category = "BP"
#' )
#'
#' # Load your custom OrgDb package outside the function.
#' library(org.Pfalciparum3D7.eg.db)
#' pfdb <- org.Pfalciparum3D7.eg.db
#' gores.org <- easytopGO(
#'   geneID = geneList, bkggset = background.gset,
#'   OrgDb = pfdb, orgdb_keytype = "GID",
#'   stats = "fisher", category = "BP"
#' )
#' # KS
#' ok.ks <- ok & !is.na(res$padj)
#' bg.ks <- res$Geneid[ok.ks]
#' geneList.ks <- stats::setNames(res$padj[ok.ks], bg.ks)
#' gores.ks <- easytopGO(
#'   geneID = geneList.ks, bkggset = bg.ks,
#'   gaf = "annotation.gaf.gz", useBiomart = FALSE, useGAF = TRUE,
#'   stats = "ks", category = "BP"
#' )
#' }
easytopGO <- function(geneID, bkggset = NULL, gaf = "",
                      useBiomart = TRUE, useGAF = FALSE,
                      mart = "protists_mart", gset = "pfalciparum_eg_gene",
                      algo = "weight01", stats = "fisher", category = "BP",
                      fdr = FALSE, correction = "BY", cutoff = 0.05,
                      go_cutoff = 0.05, return_GOdata = FALSE,
                      OrgDb = NULL, orgdb_keytype = "GID") {
  pkgs <- c("topGO", "biomaRt", "AnnotationForge")

  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_pkgs) > 0) {
    stop(
      "Required package(s) missing: ",
      paste(missing_pkgs, collapse = ", "),
      ". Please install them using BiocManager::install().",
      call. = FALSE
    )
  }


  ## 1. Validate user input before making an online BioMart request.
  stats <- match.arg(stats, c("fisher", "ks"))
  category <- match.arg(category, c("BP", "MF", "CC"))

  if (is.null(bkggset) || !length(bkggset)) {
    stop("bkggset is required: supply ALL genes eligible for the analysis")
  }
  bkggset <- unique(as.character(bkggset))
  if (anyNA(bkggset) || any(!nzchar(bkggset))) {
    stop("bkggset must contain nonmissing, nonempty gene IDs")
  }

  if (!is.numeric(geneID) || !length(geneID) || is.null(names(geneID)) ||
    anyNA(geneID) || any(!is.finite(geneID)) ||
    any(geneID < 0 | geneID > 1) || anyNA(names(geneID)) ||
    any(!nzchar(names(geneID))) || anyDuplicated(names(geneID))) {
    stop("geneID must be a uniquely named, finite numeric p-value vector in [0,1]")
  }
  if (!all(names(geneID) %in% bkggset)) {
    stop(
      "geneID contains genes absent from bkggset: ",
      paste(setdiff(names(geneID), bkggset), collapse = ", ")
    )
  }
  if (stats == "ks" && !setequal(names(geneID), bkggset)) {
    stop(
      "KS requires genuine scores for EVERY gene in bkggset; ",
      "DEG-only input is valid for Fisher, not KS"
    )
  }

  if (!is.logical(useBiomart) || length(useBiomart) != 1L ||
    is.na(useBiomart) || !is.logical(useGAF) ||
    length(useGAF) != 1L || is.na(useGAF)) {
    stop("useBiomart and useGAF must each be TRUE or FALSE")
  }
  if (!is.null(OrgDb)) {
    # An OrgDb object takes precedence over the *default* BioMart setting.
    # Reject explicitly selected competing sources instead of silently ignoring them.
    if ((!missing(useBiomart) && useBiomart) || useGAF) {
      stop("OrgDb cannot be combined with useBiomart=TRUE or useGAF=TRUE")
    }
    useBiomart <- FALSE
  } else if (identical(useBiomart, useGAF)) {
    stop("Select exactly one GO annotation source: useBiomart=TRUE or useGAF=TRUE")
  }

  # Avoid expensive online calls for unsupported test/algorithm combinations.
  algos <- if (stats == "ks") {
    c("classic", "elim", "weight01", "lea")
  } else {
    c("classic", "elim", "weight", "weight01", "lea", "parentchild")
  }
  if (!is.character(algo) || length(algo) != 1L || is.na(algo) ||
    !(algo %in% algos)) {
    stop(
      "Unsupported algorithm for ", stats, ": ", algo,
      ". Choose: ", paste(algos, collapse = ", ")
    )
  }
  if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff) ||
    !is.finite(cutoff) || cutoff <= 0 || cutoff >= 1 ||
    !is.numeric(go_cutoff) || length(go_cutoff) != 1L ||
    is.na(go_cutoff) || !is.finite(go_cutoff) ||
    go_cutoff <= 0 || go_cutoff >= 1) {
    stop("cutoff and go_cutoff must each be one number strictly between 0 and 1")
  }
  if (!is.logical(fdr) || length(fdr) != 1L || is.na(fdr) ||
    !is.logical(return_GOdata) || length(return_GOdata) != 1L ||
    is.na(return_GOdata)) {
    stop("fdr and return_GOdata must each be TRUE or FALSE")
  }
  if (!is.character(correction) || length(correction) != 1L ||
    is.na(correction) || !(correction %in% stats::p.adjust.methods)) {
    stop("correction must be one of stats::p.adjust.methods")
  }

  ## 2. Obtain GO mappings, retaining the original object name gene_2_GO.
  if (!requireNamespace("topGO", quietly = TRUE)) {
    stop("Package 'topGO' is required to run easytopGO()")
  }
  if (!is.null(OrgDb)) {
    # .usetopGO() extracts GO mappings directly from the OrgDb object.
    gene_2_GO <- NULL
  } else if (useBiomart) {
    if (!requireNamespace("biomaRt", quietly = TRUE)) {
      stop("Package 'biomaRt' is required when useBiomart = TRUE")
    }
    # Preserve the original BioMart workflow and db/go_ids variable names.
    db <- biomaRt::useMart(
      biomart = mart, host = "https://protists.ensembl.org", dataset = gset
    )
    go_ids <- biomaRt::getBM(
      attributes = c("go_id", "ensembl_gene_id"),
      filters = "ensembl_gene_id", values = bkggset, mart = db
    )
    if (!all(c("go_id", "ensembl_gene_id") %in% names(go_ids))) {
      stop("Unexpected BioMart response: missing go_id or ensembl_gene_id")
    }
    go_ids$go_id <- trimws(as.character(go_ids$go_id))
    go_ids$ensembl_gene_id <- trimws(as.character(go_ids$ensembl_gene_id))
    valid <- !is.na(go_ids$ensembl_gene_id) & nzchar(go_ids$ensembl_gene_id) &
      go_ids$ensembl_gene_id %in% bkggset & !is.na(go_ids$go_id) &
      grepl("^GO:[0-9]{7}$", go_ids$go_id)
    go_ids <- go_ids[valid, , drop = FALSE]
    gene_2_GO <- split(go_ids$go_id, go_ids$ensembl_gene_id)
  } else {
    if (!requireNamespace("readr", quietly = TRUE)) {
      stop("Package 'readr' is required when useGAF = TRUE")
    }
    # GAF col 2 = DB Object ID; col 4 = qualifier/relation; col 5 = GO ID.
    if (!is.character(gaf) || length(gaf) != 1L || is.na(gaf) || !nzchar(gaf)) {
      stop("gaf must be a valid path or URL when useGAF=TRUE")
    }
    go_ids <- readr::read_tsv(
      gaf,
      comment = "!", col_names = FALSE, trim_ws = TRUE,
      col_types = readr::cols(.default = readr::col_character()),
      progress = FALSE, show_col_types = FALSE
    )
    if (ncol(go_ids) < 5L) {
      stop("Invalid GAF: expected at least 5 tab-separated columns")
    }
    gene <- trimws(as.character(go_ids[[2L]]))
    qualifier <- trimws(as.character(go_ids[[4L]]))
    go <- trimws(as.character(go_ids[[5L]]))
    # NOT can be combined with a relation/qualifier using '|'.
    negated <- !is.na(qualifier) & grepl("(^|\\|)NOT(\\||$)", qualifier)
    # Preserve the complete non-negated GAF mapping. topGO itself intersects
    # its annotations with the universe supplied through allGenes.
    valid <- !is.na(gene) & nzchar(gene) &
      !is.na(go) & grepl("^GO:[0-9]{7}$", go) & !negated
    gene_2_GO <- split(go[valid], gene[valid])
  }

  if (is.null(OrgDb)) {
    gene_2_GO <- lapply(gene_2_GO, unique)
    if (!length(gene_2_GO) ||
      !length(intersect(names(gene_2_GO), bkggset))) {
      stop(
        "No usable GO annotations overlap bkggset. Check gene-ID namespace, ",
        "annotation file and/or BioMart dataset."
      )
    }
  }

  ## 3. Pass the same eligible-gene universe for all three annotation routes.
  goEnrichment <- .usetopGO(
    stats = stats, category = category, geneID = geneID,
    gene_2_GO = gene_2_GO, bkggset = bkggset,
    OrgDb = OrgDb, orgdb_keytype = orgdb_keytype,
    algo = algo, fdr = fdr, correction = correction,
    cutoff = cutoff, go_cutoff = go_cutoff,
    return_GOdata = return_GOdata
  )

  return(goEnrichment)
}
