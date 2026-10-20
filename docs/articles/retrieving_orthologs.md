# Retrieving Orthologs

Abstract

This is a quick and dirty tutorial on how to fetch orthologs from
OrthoMCL and InparanoiDB databases.

## Introduction

Users might require Ortholog mapping to:

- to assess phyletic profile of a set of genes and see if they have
  orthologs in other related organism or not.

- integrate datasets obtained from two different species. Often such
  datasets are scRNASeq (eg: ([Tebben et al. 2022](#ref-Tebben2022)))

- borrow annotations such as localization for hypothesis generation.

- perform co-expression analysis between two different species.

Whatever the use case maybe, we will see how we can exploit
[`getpairedOrthologs()`](https://rohit-satyam.github.io/plasmoRUtils/reference/getpairedOrthologs.md)
to fetch orthologs of multiple organisms of interest in batch manner and
integrate them.

With release 71 of VEuPathDB, only users with subscription and API key
can access the database programmatically. Use
`usethis::edit_r_environ(scope = "user")` to add VEuPathDB API using
`VEUPATHDB_API_KEY` variable. By default, all the functions look for
`VEUPATHDB_API_KEY` variable for API key. Should you store your key
using any other variable, supply it using `api_key` argument in
VEuPathDB associated functions.

``` r

# Load package and some other useful packages by using
suppressPackageStartupMessages(
  suppressWarnings({
    library(plasmoRUtils)
    library(dplyr)}))
```

Here, we will try to get all *Plasmodium falciparum* 3D7 orthologs
present in *Toxoplasma gondii* ME49, *Plasmodium berghei* ANKA and
*Plasmodium vivax* Sal1.

``` r


#?getpairedOrthologs

ids <- listOrthomcl()

numids <- c(
  ids %>% dplyr::filter(Organism == "Plasmodium vivax Sal-1") %>% dplyr::pull(ID),
  ids %>% dplyr::filter(Organism == "Toxoplasma gondii ME49") %>% dplyr::pull(ID),
  ids %>% dplyr::filter(Organism == "Plasmodium berghei ANKA") %>% dplyr::pull(ID)
)

res <- lapply(numids, function(x) {
  getpairedOrthologs(
    from = ids %>%
      dplyr::filter(Organism == "Plasmodium falciparum 3D7") %>%
      dplyr::pull(ID),
    to = x,
    db = "orthomcl",
    transform = FALSE
  )
}) %>% setNames(c("Pvivax", "Tgme49", "Pb"))

# Merge the three results
merged_df <- Reduce(
  function(x, y) merge(x, y, by = "Accession", all = TRUE),
  res
)

# Remove organism prefixes
strip_prefix <- function(x) sub("^[^|]+\\|", "", x)

clean_list <- function(s) {
  if (is.na(s)) return(NA_character_)

  toks <- trimws(strsplit(s, ",", fixed = TRUE)[[1]])
  paste(vapply(toks, strip_prefix, character(1)), collapse = ",")
}

# Clean only target-organism columns
merged_df[-1] <- lapply(
  merged_df[-1],
  function(col) vapply(col, clean_list, character(1))
)

# Rename columns
colnames(merged_df) <- c(
  "OrthoMCL Accession",
  "P vivax Sal1",
  "T gondii ME49",
  "P berghei ANKA"
)

head(merged_df)
#>   OrthoMCL Accession      P vivax Sal1     T gondii ME49    P berghei ANKA
#> 1              10181 No Ortholog Found No Ortholog Found No Ortholog Found
#> 2              10190 No Ortholog Found No Ortholog Found No Ortholog Found
#> 3              10214 No Ortholog Found No Ortholog Found No Ortholog Found
#> 4            1064799 No Ortholog Found No Ortholog Found No Ortholog Found
#> 5            1067437 No Ortholog Found No Ortholog Found No Ortholog Found
#> 6            1067464 No Ortholog Found No Ortholog Found No Ortholog Found
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
#>   Organism                  project_id Species            Species NCBI taxon I…¹
#>   <chr>                     <chr>      <chr>                               <dbl>
#> 1 Plasmodium falciparum 3D7 PlasmoDB   Plasmodium falcip…                   5833
#> 2 Plasmodium vivax Sal-1    PlasmoDB   Plasmodium vivax                     5855
#> 3 Toxoplasma gondii ME49    ToxoDB     Toxoplasma gondii                    5811
#> # ℹ abbreviated name: ¹​`Species NCBI taxon ID`
df <- lapply(1:nrow(dbs), function(x){
  plasmoRUtils::getTable(org=dbs[x,]$Organism, db=tolower(dbs[x,]$project_id),customFields = c("primary_key" ,"gene_orthomcl_name"))
})

## Retain protein coding genes

df2 <- lapply(df, function(x){
  x[!(stringr::str_detect(pattern = "N/A",string = x$`Ortholog Group`)),]
})

## Combine all the tables
merged_df2 <- Reduce(function(x, y) merge(x, y, by = "Ortholog Group", all = TRUE), df2)

merged_df2 %>% head(n = 10)
#>    Ortholog Group     Gene ID.x  Gene ID.y       Gene ID
#> 1     OG8_0000001 PF3D7_1315100 PVX_122575 TGME49_243500
#> 2     OG8_0000003 PF3D7_1123100 PVX_095165 TGME49_237890
#> 3     OG8_0000003 PF3D7_1123100 PVX_095165 TGME49_255710
#> 4     OG8_0000003 PF3D7_1123100 PVX_095165 TGME49_262540
#> 5     OG8_0000003 PF3D7_1123100 PVX_095165 TGME49_305860
#> 6     OG8_0000003 PF3D7_1123100 PVX_095165 TGME49_225490
#> 7     OG8_0000003 PF3D7_1123100 PVX_095165 TGME49_235750
#> 8     OG8_0000003 PF3D7_1123100 PVX_095165 TGME49_239420
#> 9     OG8_0000003 PF3D7_1123100 PVX_095165 TGME49_316150
#> 10    OG8_0000003 PF3D7_1123100 PVX_095165 TGME49_301440
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
#>  1 OG8_0000001     PF3D7_1315100              PVX_122575         TGME49_243500  
#>  2 OG8_0000003     PF3D7_1123100,PF3D7_07045… PVX_095165,PVX_08… TGME49_237890,…
#>  3 OG8_0000008     PF3D7_1454300,PF3D7_14413… PVX_117765,PVX_11… TGME49_291050,…
#>  4 OG8_0000011     PF3D7_0517300,PF3D7_09330… PVX_091955,PVX_09… TGME49_208970,…
#>  5 OG8_0000016     NA                         NA                 TGME49_226300  
#>  6 OG8_0000023     PF3D7_0922900              PVX_099555         TGME49_217740  
#>  7 OG8_0000045     NA                         NA                 TGME49_306450  
#>  8 OG8_0000051     NA                         NA                 TGME49_229140,…
#>  9 OG8_0000056     PF3D7_0422000              PVX_090140         TGME49_271888  
#> 10 OG8_0000060     NA                         NA                 TGME49_234570,…
```

## Fetching orthologs from InParanoiDB

InParanoiDB uses NCBI taxonomy IDs as organism identifiers. For *P.
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

## Fetching Preconfigured tables from OrthoMCL-DB and map old OG IDs to new IDs

OrthoMCL gets frequently updated and with every update come new Ortholog
IDs. When using published dataset, its highly likely that users would
want to map the old orthogroups to new orthogroups or assess which
proteins changed the orthogroups and got clubbed, which orthogroup got
split etc.
[`getPreconfiguredTableOrthomcl()`](https://rohit-satyam.github.io/plasmoRUtils/reference/getPreconfiguredTableOrthomcl.md)
can help you do this by accessing preconfigured table used by OrthoMCL
to do that.

``` r

## Fetching all OG IDs that maps to IDs of interest

## New to old IDs. Old Ids are present in second and third column
ids = c("OG7_0008348", "OG7_0003896")
getPreconfiguredTableOrthomcl(ids,customField = "previousGroups")
#> # A tibble: 5 × 3
#>   `Ortholog Group...1` `Ortholog Group...2` `Previous Ortholog Groups`
#>   <chr>                <chr>                <chr>                     
#> 1 OG8_0001453          OG8_0001453          OG7_0008348               
#> 2 OG8_0002935          OG8_0002935          OG7_0003896               
#> 3 OG8_0010463          OG8_0010463          OG7_0008348               
#> 4 OG8_0010910          OG8_0010910          OG7_0008348               
#> 5 OG8_0033410          OG8_0033410          OG7_0008348

## Old Id. New IDs are present in first column
getPreconfiguredTableOrthomcl("OG3_10277")
#> # A tibble: 1 × 3
#>   `Ortholog Group...1` `Ortholog Group...2` `Previous Ortholog Groups`
#>   <chr>                <chr>                <chr>                     
#> 1 OG8_0002935          OG8_0002935          OG3_10277
```

By changing `customField`, you can also access other preconfigured
tables. For example, keyword frequencies and their associated
orthogroups can be obtained using

``` r

getPreconfiguredTableOrthomcl(customField = "KeywordFrequency") 
#> # A tibble: 147,760 × 3
#>    `Ortholog Group` keyword                                frequency
#>    <chr>            <chr>                                      <dbl>
#>  1 OG8_0000000      kinase                                      0.9 
#>  2 OG8_0000000      protein kinase                              0.7 
#>  3 OG8_0000000      serine/threonine protein kinase             0.2 
#>  4 OG8_0000001      chromosome undetermined                     0.5 
#>  5 OG8_0000001      kinase                                      0.5 
#>  6 OG8_0000001      serine/Threonine kinase domain protein      0.25
#>  7 OG8_0000001      whole genome shotgun sequence               0.5 
#>  8 OG8_0000002      cbl interacting                             0.39
#>  9 OG8_0000002      cbl interacting protein kinase              0.21
#> 10 OG8_0000002      kinase                                      0.99
#> # ℹ 147,750 more rows
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
#> [1] dplyr_1.2.1        plasmoRUtils_1.1.3 rlang_1.3.0        readr_2.2.0       
#> [5] janitor_2.2.1      BiocStyle_2.32.1  
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6        xfun_0.62           bslib_0.12.0       
#>  [4] ggplot2_4.0.3       httr2_1.3.0         htmlwidgets_1.6.4  
#>  [7] websocket_1.4.4     processx_3.9.0      tzdb_0.5.0         
#> [10] vctrs_0.7.3         tools_4.4.1         generics_0.1.4     
#> [13] curl_8.1.0          parallel_4.4.1      easyPubMed_3.1.6   
#> [16] tibble_3.3.1        pkgconfig_2.0.3     data.table_1.18.6.1
#> [19] RColorBrewer_1.1-3  S7_0.2.2            desc_1.4.3         
#> [22] gt_1.4.0            lifecycle_1.0.5     compiler_4.4.1     
#> [25] farver_2.1.2        stringr_1.6.0       textshaping_1.0.5  
#> [28] chromote_0.5.1      codetools_0.2-20    snakecase_0.11.1   
#> [31] htmltools_0.5.9     sass_0.4.10         yaml_2.3.12        
#> [34] plotly_4.12.1       crayon_1.5.3        pillar_1.11.1      
#> [37] pkgdown_2.2.1       later_1.4.8         jquerylib_0.1.4    
#> [40] tidyr_1.3.2         cachem_1.1.0        tidyselect_1.2.1   
#> [43] rvest_1.0.5         digest_0.6.39       stringi_1.8.9      
#> [46] reshape2_1.4.5      purrr_1.2.2         bookdown_0.48      
#> [49] fastmap_1.2.0       grid_4.4.1          cli_3.6.6          
#> [52] magrittr_2.0.5      utf8_1.2.6          dichromat_2.0-1    
#> [55] withr_3.0.3         scales_1.4.0        promises_1.5.0     
#> [58] bit64_4.8.6         lubridate_1.9.5     timechange_0.4.0   
#> [61] rmarkdown_2.32      httr_1.4.9          bit_4.6.0          
#> [64] otel_0.2.0          ragg_1.5.2          hms_1.1.4          
#> [67] evaluate_1.0.5      knitr_1.52          bio3d_2.4-5        
#> [70] viridisLite_0.4.3   mgsub_2.0.0         Rcpp_1.1.2         
#> [73] glue_1.8.1          BiocManager_1.30.27 xml2_1.6.0         
#> [76] vroom_1.7.1         rstudioapi_0.19.0   jsonlite_2.0.0     
#> [79] plyr_1.8.9          R6_2.6.1            systemfonts_1.3.2  
#> [82] fs_2.1.0
```

## References

Tebben, Kieran, Aliou Dia, and David Serre. 2022. “Determination of the
Stage Composition of *Plasmodium* Infections from Bulk Gene Expression
Data.” *mSystems* 7 (4). <https://doi.org/10.1128/msystems.00258-22>.
