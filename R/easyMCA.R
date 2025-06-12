#' Fetch data from Malaria Cell Atlas
#'
#' This function fetches and loads the metadata and expression matrices from the desired data sets available on Malaria Cell Atlas (MCA) Database.
#'
#' @import tibble data.table
#' @import dplyr
#' @import rvest
#'
#' @param url A url of dataset from \code{listMCA} function.
#' @param type Type of data to be fetched. Use "exp" to fetch normalized and scaled values, use "raw" to get raw counts and use "data" to get the metadata of the dataset.
#'
#' @return df A dataframe of requested data type from MCA.
#' @export
#' @examples
#' \dontrun{
#'   url <- "https://www.malariacellatlas.org/downloads/pf-ch10x-set4-biorxiv.zip"
#'   # Use the function to read metadata, expression, or raw data
#'   metadata <- easyMCA(url, type = "data")
#'   expression <- easyMCA(url, type = "exp")
#'   raw_counts <- easyMCA(url, type = "raw")
#'
#'   ## make Seurat Object easily now
#'
#'   testmca <- Seurat::CreateSeuratObject(counts = raw_counts, meta.data = metadata, project = "MCA")
#' }
#'
easyMCA <- function(url, type = "data") {
  # Download and unzip the file
  zip_file <- basename(url)
  download.file(url, zip_file)
  unzip(zip_file)

  # Identify the target file based on the type
  all_files <- unzip(zip_file, list = TRUE)
  target_file <- all_files %>%
    dplyr::pull(Name) %>%
    .[grep(paste0("-", type), .)]

  if (length(target_file) == 0) stop("No file found matching the specified file type.")

  # Read the CSV file using fread
  data <- data.table::fread(target_file, sep = ",", header = TRUE) %>%
    tibble::column_to_rownames(var = names(.)[1])

  # Remove the zip and extracted files
  file.remove(c(zip_file, all_files$Name),showWarnings =FALSE)

  # Return the data
  return(data)
}

#styler:::style_active_file()
