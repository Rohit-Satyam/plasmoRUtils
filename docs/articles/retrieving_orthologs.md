# Retrieving Orthologs

Abstract

This is a quick and dirty tutorial on how to fetch orthologs from
OrthoMCL and InparanoiDB databases.

### Introduction

Users might require Ortholog mapping to:

- to assess phyletic profile of a set of genes and see if they have
  orthologs in other related organism or not.

- integrate datasets obtained from two different species. Often such
  datasets are scRNASeq (eg: ([Tebben et al. 2022](#ref-Tebben2022)))

- borrow annotations such as localization for hypothesis generation.

- perform co-expression analysis between two different species.

  Whatever the use case maybe, we will see how we can exploit
  [`getpairedOrthologs()`](https://rohit-satyam.github.io/plasmoRUtils/reference/getpairedOrthologs.md)
  to fetch orthologs of multiple organisms of interest in batch manner
  and integrate them.

``` r

# Load package and some other useful packages by using
suppressPackageStartupMessages(
  suppressWarnings({
    library(plasmoRUtils)
    library(dplyr)
    library(plyr)}))
```

Here, we will try to get all *Plasmodium falciparum* 3D7 (1742855)
orthologs present in *Toxoplasma gondii* ME49 (1747281), *Plasmodium
berghei* ANKA (1742951) and *Plasmodium vivax* Sal1(1744988). The
numbers mentioned alongside are vocabulary or unique IDs that OrthoMCL
uses to fetch paired orthologs.

``` r

#?getpairedOrthologs
# listOrthomcl()

res <- lapply(c(1744988,1747281,1742951), function(x){
  getpairedOrthologs(from = 1742855,
                     to=x,
                     db="orthomcl",transform = FALSE)
}) %>% setNames(c("Pvivax","Tgme49","Pb"))

merged_df <- Reduce(function(x, y) merge(x, y, by = "Accession", all = TRUE), res)
merged_df %>% head()
#>            Accession
#> 1 pfal|PF3D7_0100100
#> 2 pfal|PF3D7_0100300
#> 3 pfal|PF3D7_0102500
#> 4 pfal|PF3D7_0102600
#> 5 pfal|PF3D7_0102700
#> 6 pfal|PF3D7_0102800
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       Target ID.x
#> 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 pviv|PVX_110810
#> 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 pviv|PVX_110810
#> 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 pviv|PVX_110810
#> 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 pviv|PVX_088265
#> 5 pviv|PVX_112665,pviv|PVX_090250,pviv|PVX_097577,pviv|PVX_094305,pviv|PVX_101510,pviv|PVX_002500,pviv|PVX_090275,pviv|PVX_092995,pviv|PVX_090270,pviv|PVX_112675,pviv|PVX_096950,pviv|PVX_125730,pviv|PVX_101525,pviv|PVX_112670,pviv|PVX_083550,pviv|PVX_088825,pviv|PVX_092990,pviv|PVX_115465,pviv|PVX_090260,pviv|PVX_112690,pviv|PVX_088810,pviv|PVX_121897,pviv|PVX_090255,pviv|PVX_101515,pviv|PVX_112705,pviv|PVX_088850,pviv|PVX_112660,pviv|PVX_125728,pviv|PVX_109280,pviv|PVX_088820,pviv|PVX_097575,pviv|PVX_096995,pviv|PVX_090265,pviv|PVX_112655
#> 6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 pviv|PVX_081615
#>          Target ID.y
#> 1               <NA>
#> 2               <NA>
#> 3               <NA>
#> 4 tgon|TGME49_289050
#> 5               <NA>
#> 6               <NA>
#>                                                                                                                 Target ID
#> 1                                                                                                     pber|PBANKA_1332700
#> 2                                                                                                     pber|PBANKA_1332700
#> 3                                                                                                     pber|PBANKA_1332700
#> 4                                                                                                     pber|PBANKA_1225000
#> 5 pber|PBANKA_0623100,pber|PBANKA_1146100,pber|PBANKA_0524300,pber|PBANKA_0524400,pber|PBANKA_1146000,pber|PBANKA_0623300
#> 6                                                                                                     pber|PBANKA_0210300

## Tidying up the dataframe by removing organisms prefixes
strip_prefix <- function(x) sub("^[^|]+\\|", "", x)     # per-token
clean_list   <- function(s) {
    toks <- trimws(strsplit(s, ",", fixed = TRUE)[[1]])
    paste(vapply(toks, strip_prefix, character(1)), collapse = ",")
}

merged_df[] <- lapply(merged_df, function(col) vapply(col, clean_list, character(1)))

## Changing column names
colnames(merged_df) <- c("P falciparum 3D7","P vivax Sal1","T gondii ME49","Plasmodium berghei ANKA")

merged_df %>% head()
#>   P falciparum 3D7
#> 1    PF3D7_0100100
#> 2    PF3D7_0100300
#> 3    PF3D7_0102500
#> 4    PF3D7_0102600
#> 5    PF3D7_0102700
#> 6    PF3D7_0102800
#>                                                                                                                                                                                                                                                                                                                                                                            P vivax Sal1
#> 1                                                                                                                                                                                                                                                                                                                                                                            PVX_110810
#> 2                                                                                                                                                                                                                                                                                                                                                                            PVX_110810
#> 3                                                                                                                                                                                                                                                                                                                                                                            PVX_110810
#> 4                                                                                                                                                                                                                                                                                                                                                                            PVX_088265
#> 5 PVX_112665,PVX_090250,PVX_097577,PVX_094305,PVX_101510,PVX_002500,PVX_090275,PVX_092995,PVX_090270,PVX_112675,PVX_096950,PVX_125730,PVX_101525,PVX_112670,PVX_083550,PVX_088825,PVX_092990,PVX_115465,PVX_090260,PVX_112690,PVX_088810,PVX_121897,PVX_090255,PVX_101515,PVX_112705,PVX_088850,PVX_112660,PVX_125728,PVX_109280,PVX_088820,PVX_097575,PVX_096995,PVX_090265,PVX_112655
#> 6                                                                                                                                                                                                                                                                                                                                                                            PVX_081615
#>   T gondii ME49
#> 1            NA
#> 2            NA
#> 3            NA
#> 4 TGME49_289050
#> 5            NA
#> 6            NA
#>                                                                     Plasmodium berghei ANKA
#> 1                                                                            PBANKA_1332700
#> 2                                                                            PBANKA_1332700
#> 3                                                                            PBANKA_1332700
#> 4                                                                            PBANKA_1225000
#> 5 PBANKA_0623100,PBANKA_1146100,PBANKA_0524300,PBANKA_0524400,PBANKA_1146000,PBANKA_0623300
#> 6                                                                            PBANKA_0210300
```

Since OrthoMCL and other VEuPathDB database are not updated
simultaneously, chances are that you might be using old OrthoMCL ID. If
that’s the case above function will not be much helpful. The workaround
is to fetch the old OrthoMCL IDs for each organisms and their respective
databases and use it as an anchor to combine and collapse gene IDs from
all the organisms. The code snippet below demonstrate the same.

``` r

## List all organisms
list <- listVeupathdb(customFields = c("primary_key","project_id","species",'species_ncbi_tax_id'))

## Subset organisms I am interested in 
dbs <- list[grep("3D7|ME49$|Sal-1", list$Organism),]
dbs
#> # A tibble: 3 × 4
#>   Organism                  `VEuPathDB Project` Species   Species NCBI taxon I…¹
#>   <chr>                     <chr>               <chr>                      <dbl>
#> 1 Plasmodium vivax Sal-1    PlasmoDB            Plasmodi…                   5855
#> 2 Plasmodium falciparum 3D7 PlasmoDB            Plasmodi…                   5833
#> 3 Toxoplasma gondii ME49    ToxoDB              Toxoplas…                   5811
#> # ℹ abbreviated name: ¹​`Species NCBI taxon ID`
df <- lapply(1:nrow(dbs), function(x){
  plasmoRUtils::getTable(org=dbs[x,]$Organism, db=tolower(dbs[x,]$`VEuPathDB Project`),customFields = c("primary_key" ,"gene_orthomcl_name"))
})

## Retain protein coding genes

df2 <- lapply(df, function(x){
  x[!(stringr::str_detect(pattern = "N/A",string = x$`Ortholog Group`)),]
})

## Combine all the tables
merged_df2 <- Reduce(function(x, y) merge(x, y, by = "Ortholog Group", all = TRUE), df2)

merged_df2 %>% head(n = 10)
#>    Ortholog Group  Gene ID.x     Gene ID.y       Gene ID
#> 1      OG6_100000       <NA>          <NA> TGME49_231880
#> 2      OG6_100004 PVX_086920          <NA>          <NA>
#> 3      OG6_100006 PVX_083360 PF3D7_1349300 TGME49_237210
#> 4      OG6_100006 PVX_118355 PF3D7_1349300 TGME49_237210
#> 5      OG6_100007 PVX_001655          <NA>          <NA>
#> 6      OG6_100012 PVX_001960 PF3D7_1021900 TGME49_261590
#> 7      OG6_100012 PVX_001960 PF3D7_1021900 TGME49_308810
#> 8      OG6_100012 PVX_001960 PF3D7_1423100 TGME49_261590
#> 9      OG6_100012 PVX_001960 PF3D7_1423100 TGME49_308810
#> 10     OG6_100012 PVX_085325 PF3D7_1021900 TGME49_261590
## Combining all the gene IDs in each column so that we have 1 row per orthogroup
merged_df2 <- merged_df2 %>% 
  dplyr::group_by( `Ortholog Group`) %>% 
  dplyr::summarise(dplyr::across(dplyr::everything(), 
                                 ~paste(unique(.x), collapse = ",")), .groups = "drop")

## Changing column names
colnames(merged_df2) <- c("Orthogroup ID","P vivax Sal1","P falciparum 3D7","T gondii ME49")
merged_df2 %>% head(n = 10)
#> # A tibble: 10 × 4
#>    `Orthogroup ID` `P vivax Sal1`             `P falciparum 3D7` `T gondii ME49`
#>    <chr>           <chr>                      <chr>              <chr>          
#>  1 OG6_100000      NA                         NA                 TGME49_231880  
#>  2 OG6_100004      PVX_086920                 NA                 NA             
#>  3 OG6_100006      PVX_083360,PVX_118355      PF3D7_1349300      TGME49_237210  
#>  4 OG6_100007      PVX_001655                 NA                 NA             
#>  5 OG6_100012      PVX_001960,PVX_085325      PF3D7_1021900,PF3… TGME49_261590,…
#>  6 OG6_100025      PVX_081792,PVX_005565,PVX… PF3D7_1465800,PF3… TGME49_294550,…
#>  7 OG6_100026      PVX_096045,PVX_124085,PVX… PF3D7_1229100,PF3… TGME49_320460  
#>  8 OG6_100029      NA                         NA                 TGME49_207210  
#>  9 OG6_100034      PVX_089960                 NA                 NA             
#> 10 OG6_100045      PVX_095250                 PF3D7_0319700      NA
```

### Fetching orthologs from InParanoiDB

InParanoiDB uses NCBI taxon IDs as unique identifiers. For *P.
falciparum* (36329), *T. gondii* (5811), *P. vivax* (126793) and *P.
berghei* (5823), we will fetch pairwise orthologs by keeping *P.
falciparum* as constant species. We can recycle the code above as
follows.

``` r

#?getpairedOrthologs
# listipdb()

res <- lapply(c(126793,5811,5823), function(x){
  getpairedOrthologs(from = 36329,
                     to=x,
                     db="ipdb",transform = TRUE)
}) %>% setNames(c("Pvivax","Tgme49","Pb"))

res[[1]] %>% head()
#> # A tibble: 6 × 11
#>   `query_Group-id` query_Bitscore query_Species           query_Inparalog-scor…¹
#>              <dbl> <chr>          <chr>                   <chr>                 
#> 1                1 8359           Plasmodium falciparum … 1                     
#> 2                2 8119           Plasmodium falciparum … 1                     
#> 3                3 7587           Plasmodium falciparum … 1                     
#> 4                4 7155           Plasmodium falciparum … 1                     
#> 5                5 6891           Plasmodium falciparum … 1                     
#> 6                6 6507           Plasmodium falciparum … 1                     
#> # ℹ abbreviated name: ¹​`query_Inparalog-score`
#> # ℹ 7 more variables: `query_Protein-name` <chr>, `query_Seed-score` <chr>,
#> #   target_Bitscore <chr>, target_Species <chr>,
#> #   `target_Inparalog-score` <chr>, `target_Protein-name` <chr>,
#> #   `target_Seed-score` <chr>
```

Combining the results across InParanoiDB, is not as straightforward as
for OrthoMCL because of lack of unique Orthogroup ID. However, user can
use Query Uniprot ID as anchor, combine the results, convert Uniprot IDs
to gene IDs and then collapse rows since two or more Uniprot IDs might
correspond to same gene IDs.

### Fetching Preconfigured tables from OrthoMCL-DB and map old OG IDs to new IDs

OrthoMCL gets frequently updated and with every update come new Ortholog
IDs. When using published dataset, its highly likely that users would
want to map the old orthogroups to new orthogroups or assess which
proteins changed the orthogroups and got clubbed, which orthogroup got
split etc. `getPreconfiguredTableOrthomcl` can help you do this by
accessing preconfigured table used by OrthoMCL to do that.

``` r

## Fetching all OG IDs that maps to IDs of interest

## New to old IDs. Old Ids are present in second and third column
ids = c("OG7_0008348", "OG7_0003896")
getPreconfiguredTableOrthomcl(ids,customField = "previousGroups")
#> # A tibble: 10 × 3
#>    `Ortholog Group...1` `Ortholog Group...2` `Previous Ortholog Groups`
#>    <chr>                <chr>                <chr>                     
#>  1 OG7_0003896          OG7_0003896          OG6_100490                
#>  2 OG7_0003896          OG7_0003896          OG3_10277                 
#>  3 OG7_0003896          OG7_0003896          OG6_438394                
#>  4 OG7_0003896          OG7_0003896          OG4_24016                 
#>  5 OG7_0003896          OG7_0003896          OG5_126769                
#>  6 OG7_0003896          OG7_0003896          OG6_442722                
#>  7 OG7_0008348          OG7_0008348          OG6_469006                
#>  8 OG7_0008348          OG7_0008348          OG6_155101                
#>  9 OG7_0008348          OG7_0008348          OG6_534558                
#> 10 OG7_0008348          OG7_0008348          OG6_100809

## Old Id. New IDs are present in first column
getPreconfiguredTableOrthomcl("OG3_10277")
#> # A tibble: 5 × 3
#>   `Ortholog Group...1` `Ortholog Group...2` `Previous Ortholog Groups`
#>   <chr>                <chr>                <chr>                     
#> 1 OG7_0003895          OG7_0003895          OG3_10277                 
#> 2 OG7_0003896          OG7_0003896          OG3_10277                 
#> 3 OG7_0003899          OG7_0003899          OG3_10277                 
#> 4 OG7_0003902          OG7_0003902          OG3_10277                 
#> 5 OG7_0003903          OG7_0003903          OG3_10277
```

By changing `customField` arguments, you can also access other
preconfigured tables from the database. Eg. if you want keyword
frequencies and their associated orthogroups, this can be achieved using

``` r

getPreconfiguredTableOrthomcl(customField = "DomainFrequency")
#> # A tibble: 0 × 3
#> # ℹ 3 variables: Ortholog Group <chr>, keyword <chr>, frequency <chr>
```

## Session Info

``` r

utils::sessionInfo()
#> R version 4.4.1 (2024-06-14 ucrt)
#> Platform: x86_64-w64-mingw32/x64
#> Running under: Windows 11 x64 (build 26200)
#> 
#> Matrix products: default
#> 
#> 
#> locale:
#> [1] LC_COLLATE=English_India.utf8  LC_CTYPE=English_India.utf8   
#> [3] LC_MONETARY=English_India.utf8 LC_NUMERIC=C                  
#> [5] LC_TIME=English_India.utf8    
#> 
#> time zone: Asia/Riyadh
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] plyr_1.8.9         dplyr_1.2.1        plasmoRUtils_1.1.1 rlang_1.3.0       
#> [5] readr_2.2.0        janitor_2.2.1      BiocStyle_2.32.1  
#> 
#> loaded via a namespace (and not attached):
#>   [1] segmented_2.2-1             fs_2.1.0                   
#>   [3] ProtGenerics_1.36.0         matrixStats_1.5.0          
#>   [5] bitops_1.0-9                lubridate_1.9.5            
#>   [7] pRoloc_1.44.1               httr_1.4.8                 
#>   [9] RColorBrewer_1.1-3          doParallel_1.0.17          
#>  [11] tools_4.4.1                 MSnbase_2.30.1             
#>  [13] utf8_1.2.6                  R6_2.6.1                   
#>  [15] lazyeval_0.2.3              withr_3.0.3                
#>  [17] prettyunits_1.2.0           gridExtra_2.3.1            
#>  [19] preprocessCore_1.66.0       cli_3.6.6                  
#>  [21] Biobase_2.64.0              textshaping_1.0.5          
#>  [23] gt_1.3.0                    sass_0.4.10                
#>  [25] topGO_2.56.0                mvtnorm_1.4-1              
#>  [27] S7_0.2.2                    randomForest_4.7-1.2       
#>  [29] proxy_0.4-29                pkgdown_2.2.0              
#>  [31] Rsamtools_2.20.0            systemfonts_1.3.2          
#>  [33] txdbmaker_1.0.1             AnnotationForge_1.46.0     
#>  [35] dichromat_2.0-0.1           parallelly_1.48.0          
#>  [37] limma_3.60.6                rstudioapi_0.19.0          
#>  [39] impute_1.78.0               RSQLite_3.53.3             
#>  [41] FNN_1.1.4.1                 generics_0.1.4             
#>  [43] BiocIO_1.14.0               vroom_1.7.1                
#>  [45] gtools_3.9.5                dendextend_1.19.1          
#>  [47] GO.db_3.19.1                Matrix_1.7-5               
#>  [49] MALDIquant_1.22.3           drawProteins_1.24.0        
#>  [51] S4Vectors_0.42.1            abind_1.4-8                
#>  [53] lifecycle_1.0.5             yaml_2.3.12                
#>  [55] snakecase_0.11.1            SummarizedExperiment_1.34.0
#>  [57] recipes_1.3.3               SparseArray_1.4.8          
#>  [59] BiocFileCache_2.12.0        grid_4.4.1                 
#>  [61] blob_1.3.0                  promises_1.5.0             
#>  [63] crayon_1.5.3                PSMatch_1.8.0              
#>  [65] lattice_0.22-9              beachmat_2.20.0            
#>  [67] GenomicFeatures_1.56.0      annotate_1.82.0            
#>  [69] chromote_0.5.1              mzR_2.38.0                 
#>  [71] KEGGREST_1.44.1             pillar_1.11.1              
#>  [73] knitr_1.51                  GenomicRanges_1.56.2       
#>  [75] rjson_0.2.23                lpSolve_5.6.23             
#>  [77] future.apply_1.20.2         codetools_0.2-20           
#>  [79] mgsub_2.0.0                 glue_1.8.1                 
#>  [81] pcaMethods_1.96.0           data.table_1.18.4          
#>  [83] MultiAssayExperiment_1.30.3 vctrs_0.7.3                
#>  [85] png_0.1-9                   gtable_0.3.6               
#>  [87] kernlab_0.9-33              cachem_1.1.0               
#>  [89] gower_1.0.2                 xfun_0.59                  
#>  [91] prodlim_2026.03.11          S4Arrays_1.4.1             
#>  [93] coda_0.19-4.1               survival_3.8-6             
#>  [95] ncdf4_1.24                  timeDate_4052.112          
#>  [97] SingleCellExperiment_1.26.0 iterators_1.0.14           
#>  [99] hardhat_1.4.3               lava_1.9.2                 
#> [101] statmod_1.5.2               MLInterfaces_1.84.0        
#> [103] ipred_0.9-15                nlme_3.1-169               
#> [105] bit64_4.8.2                 progress_1.2.3             
#> [107] filelock_1.0.3              LaplacesDemon_16.1.8       
#> [109] GenomeInfoDb_1.40.1         bslib_0.11.0               
#> [111] affyio_1.74.0               irlba_2.3.7                
#> [113] rpart_4.1.27                otel_0.2.0                 
#> [115] colorspace_2.1-2            BiocGenerics_0.50.0        
#> [117] DBI_1.3.0                   nnet_7.3-20                
#> [119] tidyselect_1.2.1            processx_3.9.0             
#> [121] bit_4.6.0                   compiler_4.4.1             
#> [123] curl_7.1.0                  rvest_1.0.5                
#> [125] httr2_1.2.3                 graph_1.82.0               
#> [127] SparseM_1.84-2              xml2_1.6.0                 
#> [129] desc_1.4.3                  DelayedArray_0.30.1        
#> [131] plotly_4.12.0               bookdown_0.47              
#> [133] rtracklayer_1.64.0          scales_1.4.0               
#> [135] hexbin_1.28.5               affy_1.82.0                
#> [137] rappdirs_0.3.4              stringr_1.6.0              
#> [139] digest_0.6.39               mixtools_2.0.0.1           
#> [141] rmarkdown_2.31              XVector_0.44.0             
#> [143] htmltools_0.5.9             pkgconfig_2.0.3            
#> [145] SingleR_2.6.0               sparseMatrixStats_1.16.0   
#> [147] MatrixGenerics_1.16.0       dbplyr_2.6.0               
#> [149] fastmap_1.2.0               htmlwidgets_1.6.4          
#> [151] UCSC.utils_1.0.0            DelayedMatrixStats_1.26.0  
#> [153] farver_2.1.2                jquerylib_0.1.4            
#> [155] jsonlite_2.0.0              BiocParallel_1.38.0        
#> [157] mclust_6.1.3                mzID_1.42.0                
#> [159] ModelMetrics_1.2.2.2        BiocSingular_1.20.0        
#> [161] RCurl_1.98-1.19             magrittr_2.0.5             
#> [163] scuttle_1.14.0              GenomeInfoDbData_1.2.12    
#> [165] Rcpp_1.1.2                  viridis_0.6.5              
#> [167] MsCoreUtils_1.16.1          vsn_3.72.0                 
#> [169] pROC_1.19.0.1               stringi_1.8.7              
#> [171] zlibbioc_1.50.0             MASS_7.3-65                
#> [173] listenv_1.0.0               parallel_4.4.1             
#> [175] Biostrings_2.72.1           splines_4.4.1              
#> [177] hms_1.1.4                   igraph_2.3.3               
#> [179] QFeatures_1.14.2            reshape2_1.4.5             
#> [181] biomaRt_2.60.1              stats4_4.4.1               
#> [183] ScaledMatrix_1.12.0         XML_3.99-0.23              
#> [185] evaluate_1.0.5              BiocManager_1.30.27        
#> [187] tzdb_0.5.0                  foreach_1.5.2              
#> [189] tidyr_1.3.2                 purrr_1.2.2                
#> [191] future_1.70.0               clue_0.3-68                
#> [193] bio3d_2.4-5                 ggplot2_4.0.3              
#> [195] rsvd_1.0.5                  xtable_1.8-8               
#> [197] restfulr_0.0.17             AnnotationFilter_1.28.0    
#> [199] easyPubMed_3.1.6            e1071_1.7-17               
#> [201] later_1.4.8                 viridisLite_0.4.3          
#> [203] class_7.3-23                ragg_1.5.2                 
#> [205] tibble_3.3.1                websocket_1.4.4            
#> [207] memoise_2.0.1               AnnotationDbi_1.66.0       
#> [209] GenomicAlignments_1.40.0    IRanges_2.38.1             
#> [211] cluster_2.1.8.2             globals_0.19.1             
#> [213] timechange_0.4.0            caret_7.0-1                
#> [215] sampling_2.11
```

## References

Tebben, Kieran, Aliou Dia, and David Serre. 2022. “Determination of the
Stage Composition of *Plasmodium* Infections from Bulk Gene Expression
Data.” *mSystems* 7 (4). <https://doi.org/10.1128/msystems.00258-22>.
