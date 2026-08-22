#' Fetch Protein-protein interaction for given gene IDs from Hit Predict database
#'
#' This function searches the Hitpredict database to retrieve the Experimental Protein-Protein Interaction data.
#'
#' @import dplyr
#' @importFrom biomaRt useEnsemblGenomes getBM
#' @importFrom S4Vectors merge
#' @importFrom readr read_lines
#' @export
#'
#' @param geneID Single gene ID.
#' @param taxid Taxon ID of the organism of interest. Default: 36329. If taxon id of organism is not known set this to NULL.
#' @param uniprotToGID To convert Uniprot ID to gene ID. Set TRUE for \emph{Plasmodium} geneIDs only.
#' @param timeout Seconds after which timeout must be considered.
#' @param max_tries Maximum retries to access the database.
#' @return A data frame, containing 11 columns: "Interaction", "Interactor", "Name", "Experiments", "Category", "Method.Score", "Annotation.Score", "Interaction.Score", "Confidence", "QueryID", "ensembl_gene_id".
#' @examples
#' \dontrun{
#' test <- searchHP("Q8I1Q4")
#'
#' ## To use it for other organism, turn off uniprotToGID and provide taxid of the organism
#' test <- searchHP("BRCA1",taxid = "9606" , uniprotToGID = FALSE)
#' }
#'

searchHP <- function(geneID,
                     taxid = "36329",
                     uniprotToGID = TRUE,
                     timeout = 60,
                     max_tries = 3) {
  
  base_url <- "https://www.hitpredict.org"
  
  fetch_text <- function(url) {
    httr2::request(url) |>
      httr2::req_user_agent("plasmoRUtils searchHP") |>
      httr2::req_timeout(timeout) |>
      httr2::req_retry(
        max_tries = max_tries,
        retry_on_failure = TRUE,
        is_transient = function(resp) {
          httr2::resp_status(resp) %in% c(408, 429, 500, 502, 503, 504)
        }
      ) |>
      httr2::req_perform() |>
      httr2::resp_body_string()
  }
  
  blank_taxid <- is.null(taxid) ||
    length(taxid) == 0 ||
    all(is.na(taxid) | trimws(as.character(taxid)) == "")
  
  species <- if (blank_taxid) {
    "0"
  } else {
    utils::URLencode(as.character(taxid[[1]]), reserved = TRUE)
  }
  
  search_url <- paste0(
    base_url,
    "/proteins.php?Value=",
    utils::URLencode(as.character(geneID), reserved = TRUE),
    "&Species=",
    species
  )
  
  html_txt <- tryCatch(
    fetch_text(search_url),
    error = function(e) {
      warning("Could not open HitPredict search page for ", geneID, ": ", conditionMessage(e))
      return(NULL)
    }
  )
  
  if (is.null(html_txt)) {
    return(tibble::tibble())
  }
  
  webpage <- rvest::read_html(html_txt)
  
  htp_link <- webpage |>
    rvest::html_nodes("a") |>
    rvest::html_attr("href") |>
    stringr::str_subset("htp_int") |>
    unique()
  
  if (length(htp_link) == 0) {
    message("No interaction found in HitPredict Database")
    return(tibble::tibble())
  }
  
  read_hp_table <- function(link) {
    txt_path <- stringr::str_replace(
      link,
      "^\\./htp_int",
      "htp_int_txt"
    )
    
    txt_url <- paste0(base_url, "/", txt_path)
    
    txt <- fetch_text(txt_url)
    
    lines <- strsplit(txt, "\r\n|\n", perl = TRUE)[[1]]
    
    skip_lines_start <- 3
    skip_lines_end <- 1
    
    if (length(lines) <= skip_lines_start + skip_lines_end + 1) {
      return(tibble::tibble())
    }
    
    table_lines <- lines[
      (skip_lines_start + 1):(length(lines) - skip_lines_end)
    ]
    
    table_lines <- table_lines[nzchar(table_lines)]
    
    if (length(table_lines) < 2) {
      return(tibble::tibble())
    }
    
    utils::read.table(
      text = paste(table_lines, collapse = "\n"),
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ) |>
      tibble::as_tibble() |>
      dplyr::mutate(QueryID = geneID)
  }
  
  alldf <- purrr::map(htp_link, function(x) {
    tryCatch(
      {
        message(paste0("\033[0;32mPPI found for: ", geneID, "\033[0m\n"))
        read_hp_table(x)
      },
      error = function(e) {
        warning("Skipping HitPredict link for ", geneID, ": ", conditionMessage(e))
        tibble::tibble()
      }
    )
  })
  
  data <- dplyr::bind_rows(alldf)
  
  if (nrow(data) == 0) {
    message("No readable interaction table found in HitPredict Database")
    return(tibble::tibble())
  }
  
  if (uniprotToGID) {
    converted <- toGeneid(
      unique(data$Interactor),
      from = "uniprot",
      "ensembl",
      org = "Plasmodium falciparum 3D7",
      db = "plasmodb",
      customFields = c("primary_key", "uniprot_ids")
    )
    
    data <- S4Vectors::merge(
      data,
      converted,
      all.x = TRUE,
      all.y = FALSE,
      by.x = "Interactor",
      by.y = "UniProt ID(s)"
    )
  }
  
  return(data)
}
