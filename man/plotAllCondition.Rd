% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotAllCondition.R
\name{plotAllCondition}
\alias{plotAllCondition}
\title{Making plots similar to malaria.tools}
\usage{
plotAllCondition(geneID, returnData = FALSE, plotify = FALSE)
}
\arguments{
\item{geneID}{Single Gene ID of Plasmodium falciparum or Plasmodium berghi.}

\item{returnData}{Logical. Use true to return dataframe used for making plots.}

\item{plotify}{To make plots interactive using plotly.}
}
\value{
A plot (or data theirof) of TPM values across multiple stages of parasite.
}
\description{
This function retrieves data from malaria.tools and generates expression value plots (in TPM) similar to those produced by the website. Use this function to create publication-ready plots.
}
\examples{
\dontrun{
  #'   ## To get Plot similar to malaria.tools
  res <- plotAllCondition(geneID = "PBANKA_0100600")
}

}
