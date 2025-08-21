#' Fetch data tables from PhenoPlasm database
#'
#' This function searches the Phenotypes of the gene IDs in Phenoplasm database and enables users to fetch sub-tables such as Disruptability and Mutant phenotypes.
#'
#' @importFrom plyr ldply
#' @importFrom glue glue
#' @import rvest
#' @export
#'
#' @param geneID Character vector of Gene IDs.
#' @param org Abbreviation of the organism. Default "pf"
#' @param fetch Numeric. Use 1 to fetch the "Disruptability" table and 2 to fetch "Mutant phenotypes" table.
#' \itemize{
#' \emph{Plasmodium Species:}
#' \item pb: \emph{Plasmodium berghii}
#' \item pk: \emph{Plasmodium knowlesi}
#' \item pf: \emph{Plasmodium falciparum}
#' \item pc: \emph{Plasmodium chabaudi}
#' \item py: \emph{Plasmodium yoelii}
#' }
#' @return A data frame.
#' @examples
#' \dontrun{
#' ## get phenotype for few genes in plasmodium falciparum
#' df <- searchPhPl(geneID = c("PF3D7_0420300","PF3D7_0621000","PF3D7_0523800"), org="pf")
#' df <- searchPhPl(geneID = c("PF3D7_0420300","PF3D7_0621000","PF3D7_0523800"), org="pf", fetch=2)
#'
#' }
#'
searchPhPl <- function(geneID="",org="pf",fetch=1){
  pb <- "https://phenoplasm.org/advanced.php?text=&genes=&primespecies=P.%20berghei%20ANKA&approach1=3&includeapproach1=on&approach2=2&includeapproach2=on&approach3=1&includeapproach3=on&approach4=4&includeapproach4=on&approach5=6&includeapproach5=on&approach6=5&includeapproach6=on&approach7=7&includeapproach7=on&display=phenos&type=web"
  pf <- "https://phenoplasm.org/advanced.php?text=&genes=&primespecies=P.%20falciparum%203D7&approach1=3&includeapproach1=on&approach2=2&includeapproach2=on&approach3=1&includeapproach3=on&approach4=4&includeapproach4=on&approach5=6&includeapproach5=on&approach6=5&includeapproach6=on&approach7=7&includeapproach7=on&display=phenos&type=web"
  pc <- "https://phenoplasm.org/advanced.php?text=&genes=&primespecies=P.%20chabaudi%20chabaudi&approach1=3&includeapproach1=on&approach2=2&includeapproach2=on&approach3=1&includeapproach3=on&approach4=4&includeapproach4=on&approach5=6&includeapproach5=on&approach6=5&includeapproach6=on&approach7=7&includeapproach7=on&display=phenos&type=web"
  pk <- "https://phenoplasm.org/advanced.php?text=&genes=&primespecies=P.%20knowlesi%20strain%20H&approach1=3&includeapproach1=on&approach2=2&includeapproach2=on&approach3=1&includeapproach3=on&approach4=4&includeapproach4=on&approach5=6&includeapproach5=on&approach6=5&includeapproach6=on&approach7=7&includeapproach7=on&display=phenos&type=web"
  py <- "https://phenoplasm.org/advanced.php?text=&genes=&primespecies=P.%20yoelii%20yoelii%2017X&approach1=3&includeapproach1=on&approach2=2&includeapproach2=on&approach3=1&includeapproach3=on&approach4=4&includeapproach4=on&approach5=6&includeapproach5=on&approach6=5&includeapproach6=on&approach7=7&includeapproach7=on&display=phenos&type=web"

  ## Checking if phenoplasm have the user supplied ids
  temp <- (rvest::read_html(get(org)) %>% rvest::html_table())[[1]]

  if(all(length(geneID[!geneID %in% temp$Gene])>0 & geneID !="")){
    notfound <- paste(geneID[!geneID %in% temp$Gene],collapse = ' ')
    message(glue::glue("Warning: The following entered Gene ID(s) is/are either invalid or not available in PhenoPlasm database: {notfound} \n"))
  }

  if( all(unique(geneID != "") & fetch ==1) ) {
    ## for gene IDs found in Phenoplasm query it repeatedly and sanatize the results for Disruptability table

    result <- plyr::ldply(lapply(geneID[geneID %in% temp$Gene], function(x){
      df <- (rvest::read_html(paste0("https://phenoplasm.org/singlegene.php?gene=",x)) %>% rvest::html_table())[[1]]
      ## Removing special characters
      df$Reference <- gsub("\n\t","",df$Reference)
      df <- Filter(function(x)!all(is.na(x)), df) %>% .[!apply(is.na(.) | . == "", 1, all),]
      df$QueryGID <- x
      return(df)
    }))
    return(result)
  } else if ( all(unique(geneID != "") & fetch ==2) ) {
    ## for gene IDs found in Phenoplasm query it repeatedly and sanatize the results Mutant phenotypes

    result <- plyr::ldply(lapply(geneID[geneID %in% temp$Gene], function(x){
      df <- (rvest::read_html(paste0("https://phenoplasm.org/singlegene.php?gene=",x)) %>% rvest::html_table())[[2]]
      if(ncol(df)<5){ ## Sometime mutant table would be missing so in that case return Null
        message(glue::glue("Warning: The entered Gene ID {x} does not have Mutant phenotype information in PhenoPlasm database \n"))
        return(NULL)
      } else{
        ## Removing special characters
        df$Reference <- gsub("\n\t","",df$Reference)
        df <- Filter(function(x)!all(is.na(x)), df) %>% .[!apply(is.na(.) | . == "", 1, all),]
        df$QueryGID <- x
        return(df)
      }

    }))
    return(result)
  }
}
