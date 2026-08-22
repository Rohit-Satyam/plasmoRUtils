#' List species metadata present in OrthoMCL
#'
#' A convenience function to quickly fetch the species related vocabulary of ≥ 857 species used by [OrthoMCL](https://orthomcl.org/orthomcl/app) database. This function helps users choose numeric IDs of the organisms between which they wish to fetch the paired orthologs. See also: [getpairedOrthologs()].
#' 
#'
#' @import dplyr
#' @importFrom jsonlite fromJSON
#' 
#' @param api_key Character. VEuPathDB API key used for authentication.
#'   By default, the value is obtained from the `VEUPATHDB_API_KEY`
#'   environment variable.
#'   
#' @details
#' User should note that the numeric IDs used by OrthoMCL are not stable and changes upon every update. So an ID which was previously assigned to `Plasmodium falciparum 3D7` might get assigned to a new organism in new release. Therefore users are encouraged to run this function prior to running of [getpairedOrthologs()]
#' 
#'   
#' @export
#'
#' @return A data frame containing information about species present in OrthoMCL and their database associated ID.
#' @examples
#' \dontrun{
#' listOrthomcl()
#' }
#'
listOrthomcl <- function(
    api_key = Sys.getenv("VEUPATHDB_API_KEY")
) {
  
  if (!nzchar(api_key)) {
    stop("VEuPathDB API key is missing.", call. = FALSE)
  }
  
  url <- paste0(
    "https://orthomcl.org/orthomcl/service/",
    "record-types/sequence/searches/BySharedOrtholog"
  )
  
  resp <- httr2::request(url) %>%
    httr2::req_auth_bearer_token(api_key) %>%
    httr2::req_timeout(60) %>%
    httr2::req_options(
      connecttimeout = 30
    ) %>%
    httr2::req_retry(
      max_tries = 3,
      retry_on_failure = TRUE
    ) %>%
    httr2::req_perform()
  
  json_data <- httr2::resp_body_json(
    resp,
    simplifyVector = TRUE
  )
  
  vocab <- json_data$searchData$parameters$vocabulary
  
  df <- as.data.frame(
    do.call(rbind, vocab),
    stringsAsFactors = FALSE
  )
  
  colnames(df) <- c("ID", "Organism", "extra")
  
  df <- df %>%
    dplyr::select(ID, Organism) %>%
    dplyr::mutate(ID=as.numeric(ID)) %>%
    unique()
  
  return(df)
}
