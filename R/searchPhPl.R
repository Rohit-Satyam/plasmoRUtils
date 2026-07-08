#' Fetch data tables from PhenoPlasm database
#'
#' This function searches the Phenotypes of the gene IDs in Phenoplasm database and enables users to fetch sub-tables such as Disruptability and Mutant phenotypes.
#'
#' @importFrom plyr ldply
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

searchPhPl <- function(geneID = "", org = "pf", fetch = 1) {
  
  ## First tries normal rvest::read_html().
  ## If that fails, retries with Latin-1 conversion.
  read_phpl_tables <- function(url) {
    tryCatch(
      {
        rvest::read_html(url) %>%
          rvest::html_table()
      },
      error = function(e) {
        message("read_html/html_table failed; retrying with Latin-1 conversion: ", url)
        
        con <- url(url, open = "rb")
        on.exit(close(con), add = TRUE)
        
        raw_txt <- rawToChar(readBin(con, what = "raw", n = 1e8))
        txt <- iconv(raw_txt, from = "latin1", to = "UTF-8", sub = "")
        
        rvest::read_html(txt) %>%
          rvest::html_table()
      }
    )
  }
  
  pb <- "https://phenoplasm.org/advanced.php?text=&genes=&primespecies=P.%20berghei%20ANKA&approach1=3&includeapproach1=on&approach2=2&includeapproach2=on&approach3=1&includeapproach3=on&approach4=4&includeapproach4=on&approach5=6&includeapproach5=on&approach6=5&includeapproach6=on&approach7=7&includeapproach7=on&approach8=8&includeapproach8=on&approach9=9&includeapproach9=on&display=all&type=web"
  pf <- "https://phenoplasm.org/advanced.php?text=&genes=&primespecies=P.%20falciparum%203D7&approach1=3&includeapproach1=on&approach2=2&includeapproach2=on&approach3=1&includeapproach3=on&approach4=4&includeapproach4=on&approach5=6&includeapproach5=on&approach6=5&includeapproach6=on&approach7=7&includeapproach7=on&approach8=8&includeapproach8=on&approach9=9&includeapproach9=on&display=all&type=web"
  pc <- "https://phenoplasm.org/advanced.php?text=&genes=&primespecies=P.%20chabaudi%20chabaudi&approach1=3&includeapproach1=on&approach2=2&includeapproach2=on&approach3=1&includeapproach3=on&approach4=4&includeapproach4=on&approach5=6&includeapproach5=on&approach6=5&includeapproach6=on&approach7=7&includeapproach7=on&approach8=8&includeapproach8=on&approach9=9&includeapproach9=on&display=all&type=web"
  pk <- "https://phenoplasm.org/advanced.php?text=&genes=&primespecies=P.%20knowlesi%20strain%20H&approach1=3&includeapproach1=on&approach2=2&includeapproach2=on&approach3=1&includeapproach3=on&approach4=4&includeapproach4=on&approach5=6&includeapproach5=on&approach6=5&includeapproach6=on&approach7=7&includeapproach7=on&approach8=8&includeapproach8=on&approach9=9&includeapproach9=on&display=all&type=web"
  py <- "https://phenoplasm.org/advanced.php?text=&genes=&primespecies=P.%20yoelii%20yoelii%2017X&approach1=3&includeapproach1=on&approach2=2&includeapproach2=on&approach3=1&includeapproach3=on&approach4=4&includeapproach4=on&approach5=6&includeapproach5=on&approach6=5&includeapproach6=on&approach7=7&includeapproach7=on&approach8=8&includeapproach8=on&approach9=9&includeapproach9=on&display=all&type=web"
  
  ## Checking if PhenoPlasm has the user-supplied IDs
  temp <- read_phpl_tables(get(org))[[1]]
  
  if (all(length(geneID[!geneID %in% temp$Gene]) > 0 & geneID != "")) {
    notfound <- paste(geneID[!geneID %in% temp$Gene], collapse = " ")
    message(
      paste0(
        "Warning: The following entered Gene ID(s) is/are either invalid or not available in PhenoPlasm database: ",
        paste(notfound, collapse = " "),
        " \n"
      )
    )
  }
  
  if (all(unique(geneID != "") & fetch == 1)) {
    
    ## For gene IDs found in PhenoPlasm query repeatedly and sanitize Disruptability table
    result <- plyr::ldply(lapply(geneID[geneID %in% temp$Gene], function(x) {
      
      tables <- read_phpl_tables(
        paste0("https://phenoplasm.org/singlegene.php?gene=", x)
      )
      
      df <- tables[[1]]
      
      ## Removing special characters
      df$Reference <- gsub("\n\t", "", df$Reference)
      df <- Filter(function(x) !all(is.na(x)), df) %>%
        .[!apply(is.na(.) | . == "", 1, all), ]
      
      df$QueryGID <- x
      print(x)
      return(df)
    }))
    
    return(result)
    
  } else if (all(unique(geneID != "") & fetch == 2)) {
    
    ## For gene IDs found in PhenoPlasm query repeatedly and sanitize Mutant phenotypes
    result <- plyr::ldply(lapply(geneID[geneID %in% temp$Gene], function(x) {
      
      tables <- read_phpl_tables(
        paste0("https://phenoplasm.org/singlegene.php?gene=", x)
      )
      
      df <- tables[[2]]
      
      if (ncol(df) < 5) {
        ## Sometimes mutant table is missing, so return NULL
        message(
          paste0(
            "Warning: The entered Gene ID ",
            x,
            " does not have Mutant phenotype information in PhenoPlasm database \n"
          )
        )
        return(NULL)
      } else {
        
        ## Removing special characters
        df$Reference <- gsub("\n\t", "", df$Reference)
        df <- Filter(function(x) !all(is.na(x)), df) %>%
          .[!apply(is.na(.) | . == "", 1, all), ]
        
        df$QueryGID <- x
        return(df)
      }
    }))
    
    return(result)
  }
}
