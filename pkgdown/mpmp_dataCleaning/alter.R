library(dplyr)
library(tidyr)
library(rvest)
library(plasmoRUtils)
library(mgsub)
clickdf <- listmpmp %>% subset(.,clickableTable=="clickable")
nonclickdf <- listmpmp %>% subset(.,clickableTable=="nonclickable")

nonclick.res.50 <- lapply(nonclickdf$URL[1:50], .nonclickable)
nonclick.res.150 <- lapply(nonclickdf$URL[51:150], .nonclickable)
nonclick.res.202 <- lapply(nonclickdf$URL[151:202], .nonclickable)
nonclick.res <- c(nonclick.res.50,nonclick.res.150,nonclick.res.202)
names(nonclick.res) <- nonclickdf$Pathways.of.MPMP

saveRDS(nonclick.res,"nonclick.rds")
click.res <- lapply(clickdf$URL, .clickable)
names(click.res) <- clickdf$Pathways.of.MPMP
saveRDS(click.res,"click.rds")

## NonClickable processing steps
  non <- readRDS("../nonclick.rds")
  non <- non[lengths(non)!=0]
  ## For ease of processing make a dataframe
  non <- stack(non) ##33374
  ## Remove version from PFids
  non$values[grep("PF3D7",non$values)] <- GeneStructureTools::removeVersion(non$values[grep("PF3D7",non$values)])
  ## Remove white space
  non$values <- trimws(non$values)
  ## Some genes in MPMP are mentioned by gene symbols. Convert them to PF IDs
  subset.plasmodb <- subset(pfPlasmodbv68, pfPlasmodbv68$`Gene Name or Symbol` %in% non$values)
  non$values[match(subset.plasmodb$`Gene Name or Symbol`,non$values)] <- subset.plasmodb$`Gene ID`
  non <- non %>% unique()
  ## Now let's take care of IDs which are not PFIDs
  round1.map <- toPfid(unique(grep("PF3D7",non$values,invert = TRUE, value = TRUE)), from = "old","ensembl")
  matched_indices <- match(non$values, round1.map$`Previous ID(s)`)
  non$values <- ifelse(is.na(matched_indices), non$values, round1.map$`Gene ID`[matched_indices])


  ## Now we will see which IDs failed to map by converting pfids to description
  toPfid(unique(non$values),"ensembl")

  ## We see there are some strings that got picked up with the PFids. Let's remove them
  remove <- "MEELKTCPFTNVIPYGLLY|AEHLARIIGTEEDRVNCPF|AHGRINPFMSSPCHIQIIAR|MNIAIDEYGQPFVILR|MATSEELKQLRC|ATSEELK|MFFSSGFPFDSMGGQQAR|MDRQVINLEEILPF|MNEQIDPFYEAK|AQKNTNKKPFGNTLTNILCK|MDVNDNPFIK|VEKKIPFINFLQEK|SYATCPFNYVNINPNQKE|Localisation|PFID|Enzyme|Gene|annotation|E.C|New|Transcript|Surface|of|iRBCs|Surface|of|gametocyte|Surface|of|ookinete|Surface|anion|channel|SPF30|SPF27|Avg|fold|change|Old|Sporozoite|Surface|K.FSPFSNESCEFTYQNAK.Y|Classes|PFam|Description"
  non <- non[!grepl(remove,non$values),]

  ## Now we will see which IDs failed to map by converting pfids to description
  toPfid(unique(non$values),"ensembl")
  ## Remove special characters from front of pfids
  non$values <- gsub("^[^a-zA-Z0-9]+", "", non$values)

  ## some wrong IDs need correction
  unmapped <- read.csv("inst/unmapped.csv")  %>% subset(.,newids!="") %>% subset(. , x %in% non$values)
  matched_indices <- match(non$values, unmapped$x)
  non$values <- ifelse(is.na(matched_indices), non$values, unmapped$newids[matched_indices])

  ##Following genes failed to convert: PF3D7_0109900 PF3D7_0531700 PF3D7_0725700 PF3D7_0112400 PF3D7_1139000 PF3D7_0206400 PF3D7_0531900  PF3D7_0725500 PF3D7_1220800 Mal13P1.39

  ## Replacing the rifins, surfin, stevor and pfemp1 with gene ids
  pfids <- readxl::read_excel("inst/substitute_geneNames_with_pfids.xlsx")

  non$values[grep("RIFIN",non$values, ignore.case = TRUE)] <- paste0(pfids$rifins[!is.na(pfids$rifins)], collapse = ";")
  non$values[grep("STEVOR",non$values, ignore.case = TRUE)] <- paste0(pfids$stevor[!is.na(pfids$stevor)], collapse = ";")
  non$values[grep("PfEMP1",non$values, ignore.case = TRUE)] <- paste0(pfids$pfemps[!is.na(pfids$pfemps)], collapse = ";")
  non$values[grep("SURFIN",non$values, ignore.case = TRUE)] <- paste0(pfids$surfin[!is.na(pfids$surfin)], collapse = ";")

  non <- non %>% unique %>% tidyr::separate_rows(,values,sep = ";")
  saveRDS(non,"inst/nonclickable.rds")


## Clickable table
  click <- readRDS("../click.rds")
  click <- click[lengths(click)!=0]
  ## For ease of processing make a dataframe
  click <- plyr::ldply(click) ##33374
  ## Remove version from PFids
  click$PFID[grep("PF3D7",click$PFID)] <- GeneStructureTools::removeVersion(click$PFID[grep("PF3D7",click$PFID)])
  ## Remove white space
  click$PFID <- trimws(click$PFID)
  # ## Some genes in MPMP are mentioned by gene symbols. Convert them to PF IDs
  # subset.plasmodb <- toPfid(click$PFID,"old","ensembl")
  # click$values[match(subset.plasmodb$`Gene Name or Symbol`,click$values)] <- subset.plasmodb$`Gene ID`
  # click <- click %>% unique()
  ## Now let's take care of IDs which are not PFIDs
  round1.map <- toPfid(unique(click$PFID), from = "old","ensembl")
  matched_indices <- match(click$PFID, round1.map$`Previous ID(s)`)
  click$PFID <- ifelse(is.na(matched_indices), click$PFID, round1.map$`Gene ID`[matched_indices])


  ## Now we will see which IDs failed to map by converting pfids to description
  toPfid(unique(click$PFID),"ensembl")

  ## some wrong IDs need correction
  unmapped <- read.csv("inst/unmapped.csv")  %>% subset(.,newids!="") %>% subset(. , x %in% click$PFID)
  matched_indices <- match(click$PFID, unmapped$x)
  click$PFID <- ifelse(is.na(matched_indices), click$PFID, unmapped$newids[matched_indices])

  ## Now we will see which IDs failed to map by converting pfids to description
  toPfid(unique(click$PFID),"ensembl")


  ##Following genes failed to convert: 1396.pre-tRNA-Met-1:tRNA pf3d7_0203500 PF3D7_0112400 PF3D7_0725700 PF3D7_0112600 PF3D7_0531900 PF3D7_0725900 PFC10_API0056

  ## Replacing the rifins, surfin, stevor and pfemp1 with gene ids
  pfids <- readxl::read_excel("inst/substitute_geneNames_with_pfids.xlsx")

  click$PFID[grep("RIFIN",click$PFID, ignore.case = TRUE)] <- paste0(pfids$rifins[!is.na(pfids$rifins)], collapse = ";")
  click$PFID[grep("STEVOR",click$PFID, ignore.case = TRUE)] <- paste0(pfids$stevor[!is.na(pfids$stevor)], collapse = ";")
  click$PFID[grep("PfEMP1",click$PFID, ignore.case = TRUE)] <- paste0(pfids$pfemps[!is.na(pfids$pfemps)], collapse = ";")
  click$PFID[grep("SURFIN",click$PFID, ignore.case = TRUE)] <- paste0(pfids$surfin[!is.na(pfids$surfin)], collapse = ";")

  click <- click %>% unique %>% dplyr::select(c(PFID,.id,)) %>% tidyr::separate_rows(,PFID,sep = ";")
  saveRDS(click,"inst/clickable.rds")

  click <- readRDS("inst/clickable.rds")
  non <- readRDS("inst/nonclickable.rds")
  colnames(click) <- c("pfids","pathway")
  colnames(non) <- c("pfids","pathway")
  all <- rbind(click, non)
saveRDS(all,"inst/mpmp.28Aug2024")

temp <- S4Vectors::merge(all,listmpmp,by.x="pathway", by.y="Pathways.of.MPMP", sort=FALSE,all.x=TRUE )
use_data(mpmp.28Aug2024)

desc <- toPfid(unique(temp$pfids),"ensembl")
temp <- S4Vectors::merge(temp,desc,by.x="pfids", by.y="Gene ID", sort=FALSE,all.x=TRUE )
mpmp.28Aug2024 <- temp

## Numbering Maps
mpmp.28Aug2024 <- mpmp.28Aug2024 %>%
  group_by(pathway) %>%
  mutate(mapid = paste0("map", cur_group_id()))


temp1 <- mpmp.28Aug2024 %>%
  group_by(mapid) %>%
  dplyr::select(pfids, mapid) %>%
  summarise(PFID_list = list(pfids)) %>%
  tibble::deframe()

temp2 <- mpmp.28Aug2024 %>%
  group_by(mapid) %>%
  dplyr::select(pathway, mapid) %>% unique()

temp2 <- setNames(temp2[[1]], temp2[[2]])
temp1 <- temp1[names(temp2)]

mpmp_gsetPathFindRv2.3.1 <- list(gene_sets=temp1,descriptions=temp2)
names(mpmp_gsetPathFindRv2.3.1$descriptions)==names(mpmp_gsetPathFindRv2.3.1$gene_sets)
use_data(mpmp.28Aug2024)
use_data(mpmp_gsetPathFindRv2.3.1,overwrite = TRUE)

## for pathfindR version 2.4.1
dftemp <- mpmp.28Aug2024
dftemp$PFID <- paste0("pfa:",dftemp$pfids)
temp1 <- dftemp %>%
  group_by(mapid) %>%
  dplyr::select(PFID, mapid) %>%
  summarise(PFID_list = list(PFID)) %>%
  tibble::deframe()


temp2 <- dftemp %>%
  dplyr::group_by(mapid) %>%
  dplyr::select(pathway, mapid) %>% unique()

temp2 <- setNames(temp2[[1]], temp2[[2]])
temp1 <- temp1[names(temp2)]

mpmp_gsetPathFindRv2.4.1 <- list(gene_sets=temp1,descriptions=temp2)
names(mpmp_gsetPathFindRv2.4.1$descriptions)==names(mpmp_gsetPathFindRv2.4.1$gene_sets)
use_data(mpmp_gsetPathFindRv2.4.1,overwrite = TRUE)
saveRDS(mpmp_gsetPathFindRv2.3.1,"mpmp_gsetPathFindRv2.3.1.rds")
saveRDS(mpmp_gsetPathFindRv2.4.1,"mpmp_gsetPathFindRv2.4.1.rds")

##

write.csv(pfPlasmodbv68,"inst/pfPlasmodbv68.csv")
