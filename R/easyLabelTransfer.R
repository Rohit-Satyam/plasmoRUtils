#' Easy label transfer from reference data
#'
#' This function retrieves data from malaria.tools and generates a dataframe containing the Stage of Parasite in which the gene is highly expressed.
#'
#' @import dplyr
#' @import SingleR
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom scuttle logNormCounts
#' @export
#'
#' @param queryCounts A object containing the raw counts or lognormalised counts of query dataset. It can be in form of dcgMatrix. If the counts are normalised, set queryNormalised=TRUE.
#' @param refCounts description A object containing the raw counts or lognormalised counts of reference dataset. It can be in form of dcgMatrix. If the counts are normalised, set refNormalised=TRUE.
#' @param referenceMeta A dataframe containing the reference metadata.
#' @param refNormalised Logical. Use TRUE, if the counts are already normalised otherwise counts are log-normalised using Scuttle's logNormCounts function.
#' @param queryNormalised Logical. Use TRUE, if the counts are already normalised otherwise counts are log-normalised using Scuttle's logNormCounts function.
#' @param labelCol Column of the metadata that contains the desired labels to be transferred to the query.
#' @param isrefBulk Logical. Use TRUE, if the reference dataset is Bulk-RNASeq. In such cases "classic" approach is used to shortlist DEGs. If false, "wilcox" method will be used for scRNAseq reference.
#' @param ... Additional arguments that can be passed to SingleR based on user's needs such as de.n=30 or aggr.ref=TRUE. Refer to SingleR documentation for further details.
#'
#' @return A DFrame object containing transferred labels which can be directly used with other functions of SingleR such as plotScoreHeatmap().
#' @examples
#' \dontrun{
#' ## Fetching the URL
#' mcalist <- listMCA()
#' data("subudhi2020") ## reference dataset
#'
#' ## Using this reference set
#' url <- "https://www.malariacellatlas.org/downloads/pf-ch10x-set4-biorxiv.zip"
#'
#' raw_counts <- easyMCA(url,type = "raw")
#' rownames(raw_counts) <- gsub("-","_",rownames(raw_counts))
#' meta <- easyMCA(url,type="data")
#'
#' ## Retaining only Asexual stage cells and Lab isolates.
#' meta <- subset(meta, meta$STAGE_LR %in% c("ring","trophozoite","schizont") & DAY != "Field")
#' raw_counts <- raw_counts[,rownames(meta)]
#'
#' ## Filtering away field isolates and asexual stage cells
#' labels <- easyLabelTransfer(queryCounts = raw_counts,
#' refCounts = subudhi2020@assays@data$counts,
#' referenceMeta = subudhi2020@colData,
#' labelCol = "timetag", isrefBulk = TRUE)
#' }
#'


easyLabelTransfer <- function(queryCounts,refCounts,referenceMeta, refNormalised=FALSE,queryNormalised=FALSE,labelCol,isrefBulk=FALSE,...){
  if(refNormalised){
    ref.sce <- SingleCellExperiment::SingleCellExperiment(list(logcounts=refCounts),colData=DataFrame(referenceMeta))
  } else {
    ref.sce <- SingleCellExperiment::SingleCellExperiment(list(counts=refCounts),colData=DataFrame(referenceMeta)) %>% scuttle::logNormCounts()
  }

  if(queryNormalised){
    query.sce <- SingleCellExperiment::SingleCellExperiment(list(logcounts=queryCounts))
  } else {
    query.sce <- SingleCellExperiment::SingleCellExperiment(list(counts=queryCounts)) %>% scuttle::logNormCounts()
  }


  if(isrefBulk){
    SingleR::SingleR(test=query.sce, ref=ref.sce, labels=ref.sce[[labelCol]],...)
  } else {
    SingleR::SingleR(test=query.sce,  ref=ref.sce, labels=ref.sce[[labelCol]],de.method="wilcox",...)
  }
}



