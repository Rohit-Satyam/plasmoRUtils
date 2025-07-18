% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/easyLabelTransfer.R
\name{easyLabelTransfer}
\alias{easyLabelTransfer}
\title{Easy label transfer from reference data}
\usage{
easyLabelTransfer(
  queryCounts,
  refCounts,
  referenceMeta,
  refNormalised = FALSE,
  queryNormalised = FALSE,
  labelCol,
  isrefBulk = FALSE,
  ...
)
}
\arguments{
\item{queryCounts}{A object containing the raw counts or lognormalised counts of query dataset. It can be in form of dcgMatrix. If the counts are normalised, set queryNormalised=TRUE.}

\item{refCounts}{description A object containing the raw counts or lognormalised counts of reference dataset. It can be in form of dcgMatrix. If the counts are normalised, set refNormalised=TRUE.}

\item{referenceMeta}{A dataframe containing the reference metadata.}

\item{refNormalised}{Logical. Use TRUE, if the counts are already normalised otherwise counts are log-normalised using Scuttle's logNormCounts function.}

\item{queryNormalised}{Logical. Use TRUE, if the counts are already normalised otherwise counts are log-normalised using Scuttle's logNormCounts function.}

\item{labelCol}{Column of the metadata that contains the desired labels to be transferred to the query.}

\item{isrefBulk}{Logical. Use TRUE, if the reference dataset is Bulk-RNASeq. In such cases "classic" approach is used to shortlist DEGs. If false, "wilcox" method will be used for scRNAseq reference.}

\item{...}{Additional arguments that can be passed to SingleR based on user's needs such as de.n=30 or aggr.ref=TRUE. Refer to SingleR documentation for further details.}
}
\value{
A DFrame object containing transferred labels which can be directly used with other functions of SingleR such as plotScoreHeatmap().
}
\description{
This function retrieves data from malaria.tools and generates a dataframe containing the Stage of Parasite in which the gene is highly expressed.
}
\examples{
\dontrun{
## Fetching the URL
mcalist <- listMCA()
data("subudhi2020") ## reference dataset

## Using this reference set
url <- "https://www.malariacellatlas.org/downloads/pf-ch10x-set4-biorxiv.zip"

raw_counts <- easyMCA(url,type = "raw")
rownames(raw_counts) <- gsub("-","_",rownames(raw_counts))
meta <- easyMCA(url,type="data")

## Retaining only Asexual stage cells and Lab isolates.
meta <- subset(meta, meta$STAGE_LR \%in\% c("ring","trophozoite","schizont") & DAY != "Field")
raw_counts <- raw_counts[,rownames(meta)]

## Filtering away field isolates and asexual stage cells
labels <- easyLabelTransfer(queryCounts = raw_counts,
refCounts = subudhi2020@assays@data$counts,
referenceMeta = subudhi2020@colData,
labelCol = "timetag", isrefBulk = TRUE)
}

}
