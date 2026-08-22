#' Fetch Minor intron tables from MiDB
#'
#' This function provides ability to fetch the intron class data for 265 species in MiDB database. For more information refer to the [MiDB database](https://midb.pnb.uconn.edu/aboutus.php)
#'
#' @import dplyr
#' @import rvest
#' @importFrom readr read_csv
#' @export
#'
#' @param org Name of Organism. Can be obtained by loading midbSpecies data from the package
#' @param type Type of data to be fetched. Default: "intron"
#' @return df This function returns a dataframe of intron classification in Plasmodium species.
#'
#' @examples
#' \dontrun{
#' load("data/midbSpecies.rda")
#' ## Fetching intron data from MiDB for P. falciparum
#' df <- searchMidb(midbSpecies$`Available Species`[196])
#' }
#'

searchMidb <- function(org, type="intron"){
  url <- paste0("https://midb.pnb.uconn.edu/Download/",gsub(" ","_",org),"_",type,".csv")
  readr::read_csv(url)
}
