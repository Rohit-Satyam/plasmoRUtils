#' Quickly fetch data from MPMP database
#'
#' This function provides ability to fetch pathways data from MPMP. This data can then be modified and used for MPMP pathway enrichment analysis.
#'
#' @import dplyr
#' @import rvest
#' @import glue
#' @export
#'
#' @param url URL of the pathway of interest.
#' @return df This function returns a dataframe of containing the gene ID and Annotations fetched from MPMP database.
#' @examples
#' \dontrun{
#' df <- getMpmp("http://mpmp.huji.ac.il/maps/HNE_prot.html")
#' df <- getMpmp("http://mpmp.huji.ac.il/maps/14-3-3prot.html")
#' }
#'
getMpmp <- function(url) {
  # Attempt to use the clickable function
  tryCatch({
    df <- .clickable(url)
    if (is.null(df) || nrow(df) == 0) {
      stop("clickable returned a NULL or empty dataframe.")
    }
    return(df)
  }, error = function(e) {
    message("clickable function failed. Attempting nonclickable function.")
    # Attempt to use the nonclickable function if clickable fails
    tryCatch({
      df <- .nonclickable(url)
      if (is.null(df) || nrow(df) == 0) {
        stop("nonclickable returned a NULL or empty dataframe.")
      }
      return(df)
    }, error = function(e) {
      message("nonclickable function also failed. No data could be retrieved.")
      return(NULL)
    })
  })
}
