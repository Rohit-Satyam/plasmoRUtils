#' plasmoRUtils
#'
#' Hidden function to fetch the expandable tables from MPMP database
#'
#' @import rvest
#' @param url MPMP url.
#' @keywords internal

.clickable <- function(url) {
  df <- rvest::read_html(url) %>%
    rvest::html_elements("table.table-bordered.table-hover") %>%
    rvest::html_table()

  if (length(df) == 0) {
    message(paste0("No table found for the given URL: ", url, "\n"))
  } else {
    message(paste0("\033[0;32mfetched successfully: ", url, "\033[0m\n"))
    return(df[[1]])
  }
}

#' plasmoRUtils
#'
#' Hidden function to download static tables from MPMP database
#'
#' @import stringr rvest dplyr
#' @importFrom janitor row_to_names
#' @param url MPMP url.
#'
#' @keywords internal
.nonclickable <- function(url) {
  df <- rvest::read_html(url) %>%
    rvest::html_elements("tr") %>%
    rvest::html_text2() %>%
    as.data.frame() %>%
    dplyr::pull(1)
  if (length(df) < 2) {
    message(paste0("No table found for the given URL: ", url, "\n"))
  } else {
    pattern <- "PF3D7.*|^MAL*|PF.*|.pre-tRNA-|RNAzID|U5RNA|ORF|Surf|.*t000.*|.*m000.*|PF.*|3D7surf.*|Pf[0-9]+|[0-9]+\\.t00[0-9]+|1396.pre-trna-gly-1|1981.m00174|MaL13P1.80|pf14_0741"
    df <- df[grep("Annotation|PfID|PlasmoDB|Gene ID|NTac|Protein ID|description|Gene ID",
      df,
      ignore.case = TRUE
    ):length(df)] %>%
      stringr::str_replace_all("\t+", ":") %>%
      stringr::str_replace_all("^:+|:+$", "") %>%
      trimws() %>%
      purrr::discard(~ . == "") %>%
      # grep("PF3D7*|MAL*|PC*|*tRNA*|ORF|RNA|pf|Surf|*t000*|*m000*", ., ignore.case = TRUE, value = TRUE) %>%
      grep(pattern, ., ignore.case = TRUE, value = TRUE) %>%
      stringr::str_split(":", simplify = TRUE) %>%
      as.data.frame()

    ## getting only PFIds
    df <- df[, which(sapply(df, function(column) any(grepl(pattern, column))))]
    ## The the returned object has more than one column itneeds to be processed
    df <- if (is.data.frame(df)) {
      df[sapply(df, function(column) grepl(pattern, column))]
    } else {
      df
    }
    df <- df[!grepl("Annotation|PfID|PlasmoDB|Gene ID|NTac|Protein ID|description|Gene ID", df)]

    ## remove longer strings:
    df <- df[!(nchar(df) >= 25)]
    # colnames(df) <- c("PfID", "Annotation")
    message(paste0("\033[0;32mfetched successfully: ", url, "\033[0m\n"))
    return(unique(df))
  }
}

#' Fastest IfElse function: https://github.com/ICJIA/r-user-group/issues/11
#'
#' @keywords internal

.fast_ifelse2 <- function(test, yes, no) {
  stopifnot(identical(class(yes), class(no)))

  out <- rep(NA, length(test))
  out[test] <- yes
  out[!test] <- no
  class(out) <- class(yes)

  out
}

#' plasmoRUtils
#'
#' Hidden function called by [easytopGO()] to run topGO.
#'
#' @import stringr dplyr
#' @importFrom methods new
#'
#' @keywords internal
# Revised .usetopGO() for plasmoRUtils.
# Retains the original local names and final output columns.
# Requires topGO.

.usetopGO <- function(stats = "fisher", category = "BP", geneID, bkggset,
                      gene_2_GO = NULL,
                      algo = "weight01", fdr = FALSE,
                      correction = "BH", cutoff = 0.05,
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
  
  require(topGO)
  
  ## 1. Validate analysis options and gene-level input.
  stats <- match.arg(stats, c("fisher", "ks"))
  category <- match.arg(category, c("BP", "MF", "CC"))

  algos <- if (stats == "ks") {
    c("classic", "elim", "weight01", "lea")
  } else {
    c("classic", "elim", "weight", "weight01", "lea", "parentchild")
  }
  if (!is.character(algo) || length(algo) != 1L || is.na(algo) ||
    !(algo %in% algos)) {
    stop(
      "Unsupported algorithm for ", stats, ". Choose: ",
      paste(algos, collapse = ", ")
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
    stop("fdr and return_GOdata must be single TRUE/FALSE values")
  }
  if (!is.character(correction) || length(correction) != 1L ||
    is.na(correction) || !(correction %in% stats::p.adjust.methods)) {
    stop("correction must be one of stats::p.adjust.methods")
  }
  if (!is.numeric(geneID) || !length(geneID) ||
    is.null(names(geneID)) || anyNA(geneID) ||
    any(!is.finite(geneID)) || any(geneID < 0 | geneID > 1) ||
    anyNA(names(geneID)) || any(!nzchar(names(geneID))) ||
    anyDuplicated(names(geneID))) {
    stop("geneID must be a uniquely named, finite numeric p-value vector in [0,1]")
  }
  if (missing(bkggset) || is.null(bkggset) || !length(bkggset)) {
    stop("bkggset is required")
  }
  bkggset <- unique(as.character(bkggset))
  if (anyNA(bkggset) || any(!nzchar(bkggset))) {
    stop("bkggset must contain nonmissing, nonempty gene IDs")
  }
  if (!all(names(geneID) %in% bkggset)) {
    stop(
      "geneID contains genes absent from bkggset: ",
      paste(setdiff(names(geneID), bkggset), collapse = ", ")
    )
  }
  if (stats == "ks" && !setequal(names(geneID), bkggset)) {
    stop(
      "KS requires genuine p-values for EVERY gene in bkggset; ",
      "do not replace missing scores with 1"
    )
  }

  ## Accept a pre-loaded OrgDb object directly, including custom GID DBs.
  ## topGO already Depends on AnnotationDbi; 
  if (!is.null(OrgDb)) {
    if (!is.null(gene_2_GO)) {
      stop("Supply either OrgDb or gene_2_GO, not both")
    }
    if (!methods::is(OrgDb, "OrgDb")) {
      stop("OrgDb must be a loaded OrgDb object, not a package-name string")
    }
    if (!is.character(orgdb_keytype) || length(orgdb_keytype) != 1L ||
      is.na(orgdb_keytype) || !nzchar(orgdb_keytype)) {
      stop("orgdb_keytype must be one nonempty identifier type")
    }

    # Access the exports of topGO's mandatory dependency (no new import).
    db_keytypes <- getExportedValue("AnnotationDbi", "keytypes")
    db_columns <- getExportedValue("AnnotationDbi", "columns")
    db_keys <- getExportedValue("AnnotationDbi", "keys")
    db_select <- getExportedValue("AnnotationDbi", "select")

    available_keytypes <- db_keytypes(OrgDb)
    if (!(orgdb_keytype %in% available_keytypes)) {
      stop(
        "Unsupported orgdb_keytype: ", orgdb_keytype,
        ". Available key types: ", paste(available_keytypes, collapse = ", ")
      )
    }
    if (!("GO" %in% db_columns(OrgDb))) {
      stop(
        "This OrgDb has no GO column. A GFF-derived database must have ",
        "GO associations added (not only genome features) for enrichment."
      )
    }

    # Preserve the original universe; only annotation mapping is restricted.
    matched_keys <- intersect(bkggset, db_keys(OrgDb, keytype = orgdb_keytype))
    if (!length(matched_keys)) {
      stop("No bkggset IDs match the OrgDb key type '", orgdb_keytype, "'")
    }
    go_ids <- suppressMessages(db_select(
      OrgDb,
      keys = matched_keys, columns = "GO", keytype = orgdb_keytype
    ))
    if (!all(c(orgdb_keytype, "GO") %in% names(go_ids))) {
      stop("Unexpected OrgDb response: missing ID or GO column")
    }
    gene <- as.character(go_ids[[orgdb_keytype]])
    go <- trimws(as.character(go_ids[["GO"]]))
    valid <- !is.na(gene) & nzchar(gene) &
      !is.na(go) & grepl("^GO:[0-9]{7}$", go)
    if (!any(valid)) {
      stop("No valid GO annotations for matching background genes in OrgDb")
    }
    gene_2_GO <- split(go[valid], gene[valid])
  }

  if (!is.list(gene_2_GO) || is.null(names(gene_2_GO)) ||
    anyNA(names(gene_2_GO)) || any(!nzchar(names(gene_2_GO))) ||
    anyDuplicated(names(gene_2_GO))) {
    stop("gene_2_GO must be a named list with unique, nonempty gene IDs")
  }

  ## Retain full annotation mapping; topGO restricts it to the gene universe.
  gene_2_GO <- lapply(gene_2_GO, function(go) {
    go <- as.character(go)
    unique(go[!is.na(go) & grepl("^GO:[0-9]{7}$", go)])
  })
  gene_2_GO <- Filter(length, gene_2_GO)
  if (!length(intersect(names(gene_2_GO), bkggset))) {
    stop("No GO annotations overlap bkggset. Check gene-ID namespace.")
  }

  ## 2. Construct topGOdata; Fisher and KS require different inputs.
  if (stats == "fisher") {
    # Binary DEG selection over ALL background genes.
    selected <- names(geneID)[geneID < cutoff]
    all_genes <- factor(as.integer(bkggset %in% selected), levels = c(0, 1))
    names(all_genes) <- bkggset

    GOdata <- methods::new(
      "topGOdata",
      ontology = category, allGenes = all_genes,
      annot = topGO::annFUN.gene2GO, gene2GO = gene_2_GO
    )
  } else {
    # KS must retain the real p-value ordering of the complete background.
    all_genes <- geneID[match(bkggset, names(geneID))]
    names(all_genes) <- bkggset

    GOdata <- methods::new(
      "topGOdata",
      ontology = category, allGenes = all_genes,
      annot = topGO::annFUN.gene2GO, gene2GO = gene_2_GO,
      geneSel = function(allScore) !is.na(allScore) & allScore < cutoff
    )
  }

  if (return_GOdata) {
    return(GOdata)
  }
  if (!length(topGO::genes(GOdata))) {
    stop("No feasible genes remain in the selected GO ontology")
  }

  ## 3. Test enrichment; use numeric scores, not formatted GenTable strings.
  results <- topGO::runTest(GOdata, algorithm = algo, statistic = stats)
  results_p <- topGO::score(results)
  if (!length(results_p)) {
    stop("No GO terms passed the annotation/node-size checks")
  }

  # Keep the original internal result-column name: 'stats'.
  goEnrichment <- topGO::GenTable(
    GOdata,
    stats = results, orderBy = "stats",
    topNodes = length(results_p)
  )

  # GenTable can print tiny p-values as '< 1e-30'. Restore exact numeric
  # scores by GO identifier before adjustment or filtering.
  goEnrichment[["stats"]] <- unname(results_p[
    match(as.character(goEnrichment$GO.ID), names(results_p))
  ])
  goEnrichment$Term <- as.character(goEnrichment$Term)

  ## 4. Adjust ALL tested GO terms before reporting significant rows.
  if (fdr) {
    results_padj <- stats::p.adjust(results_p, method = correction)
    goEnrichment[["padj"]] <- unname(results_padj[
      match(as.character(goEnrichment$GO.ID), names(results_padj))
    ])
    goEnrichment <- goEnrichment[
      !is.na(goEnrichment[["padj"]]) & goEnrichment[["padj"]] < go_cutoff, ,
      drop = FALSE
    ]
    goEnrichment <- goEnrichment[
      order(goEnrichment[["padj"]]), ,
      drop = FALSE
    ]
  } else {
    goEnrichment <- goEnrichment[
      !is.na(goEnrichment[["stats"]]) & goEnrichment[["stats"]] < go_cutoff, ,
      drop = FALSE
    ]
    goEnrichment <- goEnrichment[
      order(goEnrichment[["stats"]]), ,
      drop = FALSE
    ]
  }
  rownames(goEnrichment) <- NULL

  ## 5. Restore complete descriptions where topGO has truncated terms.
  
  .goIdToTerm <- function(go_ids) {
  go_ids <- as.character(go_ids)
  terms <- rep(NA_character_, length(go_ids))
  names(terms) <- go_ids
  valid <- !is.na(go_ids) & nzchar(go_ids)
  if (!any(valid)) return(terms)
  annotation <- AnnotationDbi::select(
    GO.db::GO.db,
    keys = unique(go_ids[valid]),
    columns = "TERM",
    keytype = "GOID"
  )

  terms[valid] <- annotation$TERM[
    match(go_ids[valid], annotation$GOID)
  ]

  return(terms)
}

  if (nrow(goEnrichment)) {
    incomplete <- which(!is.na(goEnrichment$Term) &
      grepl("\\.\\.\\.$", trimws(goEnrichment$Term)))
    if (length(incomplete)) {
      if (requireNamespace("AnnotationDbi", quietly = TRUE)) {
        full_terms <- .goIdToTerm(
          as.character(goEnrichment$GO.ID[incomplete])
        )
        valid <- which(!is.na(full_terms) & nzchar(trimws(full_terms)))
        if (length(valid)) {
          goEnrichment$Term[incomplete[valid]] <- unname(full_terms[valid])
        }
      } else {
        warning("AnnotationDbi is not installed; install it.",
          call. = FALSE
        )
      }
    }
  }

  ## Preserve the original reversed Term-factor levels for plotting.
  goEnrichment$Term <- factor(
    goEnrichment$Term,
    levels = rev(unique(goEnrichment$Term))
  )

  ## 6. Attach GO-member gene identifiers and selected GO-member genes.
  if (nrow(goEnrichment)) {
    get.genes <- topGO::genesInTerm(GOdata, as.character(goEnrichment$GO.ID))
    significant <- topGO::sigGenes(GOdata)

    goEnrichment$associated_genes <- vapply(
      as.character(goEnrichment$GO.ID),
      function(id) paste(unique(get.genes[[id]]), collapse = ","),
      character(1)
    )
    goEnrichment$significant_genes <- vapply(
      as.character(goEnrichment$GO.ID),
      function(id) {
        paste(intersect(unique(get.genes[[id]]), significant),
          collapse = ","
        )
      },
      character(1)
    )
  } else {
    # Return a correctly typed zero-row data frame, not invisible NULL.
    goEnrichment$associated_genes <- character(0)
    goEnrichment$significant_genes <- character(0)
  }

  ## 7. Restore original public output name: fisher or ks.
  colnames(goEnrichment)[colnames(goEnrichment) == "stats"] <- stats

  attr(goEnrichment, "gene_universe") <- list(
    available = length(bkggset),
    feasible = length(topGO::genes(GOdata)),
    selected_feasible = length(topGO::sigGenes(GOdata))
  )

  return(goEnrichment)
}

cleanFun <- function(htmlString) {
  return(gsub("<.*?>", "", htmlString))
}

convert_last_letter <- function(str) {
  sub("([A-Za-z])$", "\\L\\1", str, perl = TRUE)
}

utils::globalVariables(".")

#' plasmoRUtils
#'
#' Hidden function for searchRS
#'
#' @import httr2 jsonlite
#'
#' @keywords internal
.rs_api_search_page <- function(query,
                                offset = 0,
                                article_type = "Research Article",
                                status = "all") {
  resp <- httr2::request("https://www.researchsquare.com/api/search") |>
    httr2::req_url_query(
      articleType = article_type,
      offset = offset,
      status = status,
      unified = query
    ) |>
    httr2::req_user_agent("plasmoRUtils searchRS") |>
    httr2::req_timeout(60) |>
    httr2::req_retry(
      max_tries = 3,
      retry_on_failure = TRUE,
      is_transient = function(resp) {
        httr2::resp_status(resp) %in% c(408, 429, 500, 502, 503, 504)
      }
    ) |>
    httr2::req_perform()

  txt <- httr2::resp_body_string(resp)
  out <- jsonlite::fromJSON(txt, simplifyVector = TRUE)
  out$result
}

#' plasmoRUtils
#'
#' Hidden function for searchRS
#'
#'
#' @keywords internal
.escape_regex <- function(x) {
  vapply(strsplit(x, "", fixed = TRUE), function(chars) {
    meta <- chars %in% c(
      "\\", ".", "|", "(", ")", "[", "]", "{", "}",
      "^", "$", "*", "+", "?"
    )

    paste0(ifelse(meta, paste0("\\", chars), chars), collapse = "")
  }, character(1))
}


#' plasmoRUtils
#'
#' Hidden function for searchRS
#'
#' @import purrr
#'
#' @keywords internal
.make_flexible_regex <- function(terms) {
  terms <- unique(terms[!is.na(terms) & nzchar(terms)])

  if (length(terms) == 0) {
    return("(?!)")
  }

  term_regex <- purrr::map_chr(terms, function(x) {
    parts <- unlist(strsplit(x, "[-_[:space:]]+"))
    parts <- .escape_regex(parts)

    paste(parts, collapse = "[-_[:space:]]*")
  })

  paste0(
    "(?i)(?<![A-Za-z0-9])(?:",
    paste(term_regex, collapse = "|"),
    ")(?![A-Za-z0-9])"
  )
}

#' plasmoRUtils
#'
#' Hidden function for searchRS
#'
#' @import stringr httr2 rvest
#'
#' @keywords internal
.rs_article_text <- function(url_path) {
  url <- ifelse(
    stringr::str_starts(url_path, "https?://"),
    url_path,
    paste0("https://www.researchsquare.com", url_path)
  )

  resp <- httr2::request(url) |>
    httr2::req_user_agent("plasmoRUtils searchRS") |>
    httr2::req_timeout(60) |>
    httr2::req_retry(
      max_tries = 3,
      retry_on_failure = TRUE,
      is_transient = function(resp) {
        httr2::resp_status(resp) %in% c(408, 429, 500, 502, 503, 504)
      }
    ) |>
    httr2::req_perform()

  txt <- httr2::resp_body_string(resp)

  html <- rvest::read_html(txt)
  body <- rvest::html_element(html, "body")

  rvest::html_text2(body)
}

#' plasmoRUtils
#'
#' Hidden function to get tables from VEuPathDB databases with authentication
#'
#' @import httr2 readr
#'
#' @keywords internal
#'
.read_authenticated_tsv <- function(url, api_key) {
  if (!is.character(api_key) || length(api_key) != 1L ||
      is.na(api_key) || !nzchar(trimws(api_key))) {
    stop("VEuPathDB API key is missing. Set VEUPATHDB_API_KEY.", call. = FALSE)
  }
  resp <- httr2::request(url) %>%
    httr2::req_auth_bearer_token(api_key) %>%
    httr2::req_headers(
      Accept = "text/plain"
    ) %>%
    httr2::req_perform()
  
  txt <- httr2::resp_body_string(resp)
  
  readr::read_tsv(
    I(txt),
    show_col_types = FALSE,
    progress = FALSE
  )
}
