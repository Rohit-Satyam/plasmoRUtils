#' List species metadata present in OrthoMCL
#'
#' A convenience function to quickly fetch the species related vocabulary of 713 species used by [OrthoMCL](https://orthomcl.org/orthomcl/app) database. This function helps users choose IDs of the organisms between which they wish to fetch the paired orthologs. See also: \code{getpairedOrthologs()}.
#'
#' @import dplyr
#' @importFrom jsonlite fromJSON
#' @export
#'
#' @return A data frame containing information about species present in InParanoiDB9 and their taxon ID.
#' @examples
#' \dontrun{
#' listOrthomcl()
#' }
#'
listOrthomcl <- function(){
  url <- "https://orthomcl.org/orthomcl/service/record-types/sequence/searches/BySharedOrtholog"
  json_data <- jsonlite::fromJSON(url)
  
  # extract vocabulary (it's nested)
  vocab <- json_data$searchData$parameters$vocabulary
  
  # convert to dataframe
  df <- as.data.frame(do.call(rbind, vocab), stringsAsFactors = FALSE)
  colnames(df) <- c("ID", "Organism", "extra")
  
  # keep only ID and Organism
  df <- df %>% 
    dplyr::select(ID, Organism) %>% unique()
  return(df)
}
