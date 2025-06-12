utils::globalVariables(c("Description and links", "Gene ID", "Gene Name or Symbol", "GeneID", "Previous ID(s)","interactorA","interactorB","Product Description","Title","uniprotsptrembl","malariatools","gene","condition","group","accession","pfPlasmodbv68"))

sanitize <- function(string) {
  return(gsub("<.*?>", "", string))
}
