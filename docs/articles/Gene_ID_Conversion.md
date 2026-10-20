# Accessing component databases of VEuPathDB

Abstract

This section covers the capability of `plasmoRUtils` to access various
VEuPathDb database and carry out ID conversion tasks, fetch
preconfigured data tables and few other tasks.

## Introduction

The *plasmoRUtils* package streamlines access to VEuPathDB’s family of
12 specialized databases such as ToxoDB, PlasmoDB, and PiroplasmaDB
([Amos et al. 2022](#ref-amos2022veupathdb)). It provides direct data
retrieval capabilities through VEuPathDB’s RESTful API, enabling
seamless integration of biological data into R workflows. The package
supports downloading both standard and customized data tables, making it
particularly valuable for researchers needing to combine data from
multiple sources for downstream analysis.

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

## Gene ID Conversion

A common challenge in bioinformatics involves mapping between different
identifier systems across databases. For apicomplexan research, this
might include converting between UniProt IDs, legacy gene identifiers,
and current Ensembl gene IDs. The
[`toGeneid()`](https://rohit-satyam.github.io/plasmoRUtils/reference/toGeneid.md)
function addresses this need by retrieving up-to-date annotations from
VEuPathDB databases, mapping supported source ID types to gene IDs or
vice-versa.

### Retrieving Gene Annotations and Alternative IDs

The
[`toGeneid()`](https://rohit-satyam.github.io/plasmoRUtils/reference/toGeneid.md)
function enables annotation retrieval when provided with
gene/proteins/entrez IDs. For custom field selection, use
[`getTableVeupathdb()`](https://rohit-satyam.github.io/plasmoRUtils/reference/getTableVeupathdb.md)
directly.

``` r

## Get annotations for list of geneIDs for PF3D7
getTableVeupathdb(c("PF3D7_0420300", "PF3D7_0621000"))
#> # A tibble: 2 × 11
#>   `Gene ID`    source_id `Input Order` Organism `Product Description` `Input ID`
#>   <chr>        <chr>     <chr>         <chr>    <chr>                 <chr>     
#> 1 PF3D7_04203… PF3D7_04… 1             Plasmod… AP2 domain transcrip… PF3D7_042…
#> 2 PF3D7_06210… PF3D7_06… 2             Plasmod… RNA polymerase subun… PF3D7_062…
#> # ℹ 5 more variables: `UniProt ID(s)` <chr>, `Entrez Gene ID` <chr>,
#> #   `Previous ID(s)` <chr>, `Gene Name or Symbol` <chr>, `Gene Strand` <chr>

## give miscellaneous IDs and choose custom fields
getTableVeupathdb(ids = c("PF3D7_1372300", "PBANKA_1355600","A0A143ZZR8","W7JSP8" ),
                  customFields=c("primary_key","predicted_go_component","annotated_go_function"))
#> # A tibble: 3 × 6
#>   `Gene ID`        `Input ID`      `Input Order` Organism Computed GO Componen…¹
#>   <chr>            <chr>           <chr>         <chr>    <chr>                 
#> 1 PBANKA_1355600   PBANKA_1355600  2             Plasmod… N/A                   
#> 2 PF3D7_1372300    A0A143ZZR8, PF… 1             Plasmod… N/A                   
#> 3 PfNF54_110037600 W7JSP8          4             Plasmod… membrane              
#> # ℹ abbreviated name: ¹​`Computed GO Components`
#> # ℹ 1 more variable: `Curated GO Functions` <chr>
```

Users can now convert IDs from one type to another such as gene IDs,
protein IDs or entrez ids using
[`toGeneid()`](https://rohit-satyam.github.io/plasmoRUtils/reference/toGeneid.md)
and defining database is no more required.

``` r

## Convert uniprot IDs back to gene IDs. It will also provide Product description and Gene Symbol

toGeneid(inputid = c("Q8I1N6","C6KT48"),from="uniprot",to="ensembl" )
#> # A tibble: 2 × 4
#>   `Input ID` `Output ID`   `Gene ID`     Organism                 
#>   <chr>      <chr>         <chr>         <chr>                    
#> 1 Q8I1N6     PF3D7_0420300 PF3D7_0420300 Plasmodium falciparum 3D7
#> 2 C6KT48     PF3D7_0621000 PF3D7_0621000 Plasmodium falciparum 3D7
toGeneid(inputid = c("PF3D7_0420300", "PF3D7_0621000"),
         from="ensembl", to= "uniprot" )
#> # A tibble: 2 × 4
#>   `Input ID`    `Output ID` `Gene ID`     Organism                 
#>   <chr>         <chr>       <chr>         <chr>                    
#> 1 PF3D7_0420300 Q8I1N6      PF3D7_0420300 Plasmodium falciparum 3D7
#> 2 PF3D7_0621000 C6KT48      PF3D7_0621000 Plasmodium falciparum 3D7

## OR Request selected annotation fields directly from getTable(). 
# Retrieve a custom table directly instead of passing a reduced field list to toGeneid().
getTable(org = "Toxoplasma gondii ME49", db = "toxodb",
         customFields = c("primary_key", "predicted_go_component",
                          "annotated_go_function")) %>%
  dplyr::filter(`Gene ID` %in% c("TGME49_304740", "TGME49_208030"))
#> # A tibble: 2 × 3
#>   `Gene ID`     `Computed GO Components` `Curated GO Functions`
#>   <chr>         <chr>                    <chr>                 
#> 1 TGME49_208030 extracellular region     N/A                   
#> 2 TGME49_304740 N/A                      N/A
```

> getTable() function will be phased out in future releases and will be
> replaced by getTableVeupathdb() in future.

This functionality proves particularly valuable for enhancing
differential expression analysis results with comprehensive annotations,
enabling complete workflow automation without leaving command-line
interfaces on HPC systems.

With
[`toGeneid()`](https://rohit-satyam.github.io/plasmoRUtils/reference/toGeneid.md)
function, one-to-many relationships are retained in long format by
default, allowing all possible mappings to remain visible. The decision
about how to handle such cases is therefore left to the user, rather
than being imposed by the function.

For instance, a UniProt accession can sometimes be associated with
multiple gene IDs. In such cases, retain the one-to-many relationship
for inspection. If a downstream step requires a single row per UniProt
accession, collapse the distinct gene IDs into a delimited list rather
than arbitrarily discarding alternatives. Do not select an arbitrary
mapping without additional biological evidence. Because there is no
universally appropriate rule for resolving this type of ambiguity, we
avoid making this choice for the end user.

``` r

toGeneid(c("Q8I0P6" ,"A0A2P1JI60"), from = "uniprot",to = "ensembl")
#> # A tibble: 36 × 4
#>    `Input ID` `Output ID`       `Gene ID`         Organism                      
#>    <chr>      <chr>             <chr>             <chr>                         
#>  1 Q8I0P6     PF3D7_1357000     PF3D7_1357000     Plasmodium falciparum 3D7     
#>  2 Q8I0P6     PF3D7_1357100     PF3D7_1357100     Plasmodium falciparum 3D7     
#>  3 A0A2P1JI60 PF3D7_1357100     PF3D7_1357100     Plasmodium falciparum 3D7     
#>  4 A0A2P1JI60 Pf7G8-2_000459500 Pf7G8-2_000459500 Plasmodium falciparum 7G8 2019
#>  5 A0A2P1JI60 Pf7G8-2_000459600 Pf7G8-2_000459600 Plasmodium falciparum 7G8 2019
#>  6 A0A2P1JI60 Pf7G8_130061700   Pf7G8_130061700   Plasmodium falciparum 7G8     
#>  7 A0A2P1JI60 PfCD01_130062700  PfCD01_130062700  Plasmodium falciparum CD01    
#>  8 A0A2P1JI60 PfCD01_130062800  PfCD01_130062800  Plasmodium falciparum CD01    
#>  9 A0A2P1JI60 PfDd2_130063000   PfDd2_130063000   Plasmodium falciparum Dd2     
#> 10 A0A2P1JI60 PfGA01_130063200  PfGA01_130063200  Plasmodium falciparum GA01    
#> # ℹ 26 more rows
```

### Accessing Preconfigured Data Tables from VEuPathDB’s component sites

Since VEuPathDB API documentation specifically encourages to use
specific organism database as quoted below, we developed
[`getTable()`](https://rohit-satyam.github.io/plasmoRUtils/reference/getTable.md)
function to fetch fields of interests from each database separately.

> There are 12 component sites and one portal: VEuPathDB.org. The
> component sites are: AmoebaDB, CryptoDB, FungiDB, GiardiaDB, HostDB,
> MicrosporidiaDB, PiroplasmaDB, PlasmoDB, ToxoDB, TrichDB, TriTrypDB
> and VectorBase. For most record types (all but dataset and organism),
> when running a search, the portal reaches out to component sites to
> get the search results. That means it will be faster to use a
> component site directly when you can.

Some frequently required fields have been provided in the help section
of
[`getTable()`](https://rohit-satyam.github.io/plasmoRUtils/reference/getTable.md).

``` r

## To fetch table for all the genes present in an organism

getTable(org="Plasmodium falciparum 3D7", db="plasmodb",
         api_key = Sys.getenv("VEUPATHDB_API_KEY")) %>% head()
#> # A tibble: 6 × 10
#>   `Gene ID`     `Product Description`        `Gene Strand` `Gene Name or Symbol`
#>   <chr>         <chr>                        <chr>         <chr>                
#> 1 PF3D7_0100100 erythrocyte membrane protei… forward       VAR                  
#> 2 PF3D7_0100200 rifin                        reverse       RIF                  
#> 3 PF3D7_0100300 erythrocyte membrane protei… reverse       VAR                  
#> 4 PF3D7_0100400 rifin                        forward       RIF                  
#> 5 PF3D7_0100500 erythrocyte membrane protei… reverse       N/A                  
#> 6 PF3D7_0100600 rifin                        reverse       RIF                  
#> # ℹ 6 more variables: `Previous ID(s)` <chr>, `Entrez Gene ID` <chr>,
#> #   `UniProt ID(s)` <chr>, `Protein Length` <chr>, `# TM Domains` <chr>,
#> #   `SignalP 6.0 Signal Probability` <chr>

## User can also provide custom fields. For example we wish to download the Transcript abundance in Transcripts per Million (TPM) [Dataset: Transcriptome in severe vs uncomplicated malaria]

getTable(org="Plasmodium falciparum 3D7", db="plasmodb", customFields = c("primary_key","pan_24905_ns_2917")) %>% head()
#> # A tibble: 6 × 2
#>   `Gene ID`     `antisense - SFM-3_severe`
#>   <chr>                              <dbl>
#> 1 PF3D7_0100100                       0.16
#> 2 PF3D7_0100200                       0   
#> 3 PF3D7_0100300                       0.13
#> 4 PF3D7_0100400                       0   
#> 5 PF3D7_0100500                       0   
#> 6 PF3D7_0100600                       0

## The attributes supported by customFields can be inspected using
attrs <- getTableAttributes(
   db = "plasmodb",
  record_type = "transcript", api_key=Sys.getenv("VEUPATHDB_API_KEY")
)

head(attrs)
#>             name    displayName columnDataType isDisplayable isSortable
#> 1    primary_key        Gene ID         STRING          TRUE       TRUE
#> 2       overview       Overview         STRING          TRUE       TRUE
#> 3    gbrowseLink    gbrowseLink         STRING          TRUE       TRUE
#> 4        JBrowse        JBrowse          OTHER          TRUE      FALSE
#> 5 gene_source_id gene_source_id         STRING         FALSE       TRUE
#> 6      source_id      source_id         STRING         FALSE       TRUE
#>   isInReport                                 help
#> 1       TRUE                                 <NA>
#> 2      FALSE                                 <NA>
#> 3      FALSE                                 <NA>
#> 4      FALSE View this Gene in the Genome Browser
#> 5       TRUE                                 <NA>
#> 6       TRUE                                 <NA>
```

For more information on fields that can be supplied to
[`getTable()`](https://rohit-satyam.github.io/plasmoRUtils/reference/getTable.md),
use the following steps:

1.  Go to database of your interest (Say “PlasmoDB”)

2.  Click on `Annotation, curation and identifiers` tab on your left and
    select `List of IDs`

3.  Scroll down and click on
    `Build a Web Services URL from this Search >>` hyperlink.

4.  Under the section `Choose Columns:` choose fields of your interest.
    Most of these fields are included in help section of
    [`getTable()`](https://rohit-satyam.github.io/plasmoRUtils/reference/getTable.md).
    However, fields specific to a particular database such as dataset
    related fields (starts with “pan\_”) are excluded.

5.  Once you select fields, they are updated in **POST** section of the
    webpage query builder.

The VEuPathDB component databases also provide preconfigured tables that
cannot be retrieved through
[`getTable()`](https://rohit-satyam.github.io/plasmoRUtils/reference/getTable.md).
To enable users to fetch such tables, we wrote another function called
[`getPreconfiguredTable()`](https://rohit-satyam.github.io/plasmoRUtils/reference/getPreconfiguredTable.md)
the usage of which has been shown below.

``` r

## To get list of preconfigured tables

preconf <- getPreconfiguredTable(listtables = TRUE)

head(preconf)
#> [1] "3dPreds"                                                           
#> [2] "afumA1163_TF_KO_Image_NAFeaturePhenotypeImage_RSRC_PhenotypeImages"
#> [3] "afumAf293_TF_KO_Image_NAFeaturePhenotypeImage_RSRC_PhenotypeImages"
#> [4] "Alias"                                                             
#> [5] "AllProducts"                                                       
#> [6] "AlphaFoldLinkouts"

## Fetch pathway table for all the genes from MPMP database

getPreconfiguredTable(org = "Plasmodium falciparum 3D7",db = "plasmodb",customField = "MetabolicPathwaysMPMP") %>% head()
#> # A tibble: 6 × 4
#>   `Gene ID`     pathway_id          Pathway                             Activity
#>   <chr>         <chr>               <chr>                               <chr>   
#> 1 PF3D7_0100100 cytoadherencescheme Interactions between modified host… erythro…
#> 2 PF3D7_0100100 exportome           Exportome - compiled from various … erythro…
#> 3 PF3D7_0100100 gene_lumef          Gene expression affected by lumefa… erythro…
#> 4 PF3D7_0100100 lys_met             Peptides with confirmed methylated… erythro…
#> 5 PF3D7_0100100 Par_RBC             Protein-Protein Interactions betwe… erythro…
#> 6 PF3D7_0100100 PfEMP1              PfEMP1 domain architectures         erythro…
```

Please note that the MPMP pathway version provided by PlasmoDB is
[outdated](https://plasmodb.org/plasmo/app/record/dataset/DS_1e177b728b)
(`03-2019`). Some pathways have been revised or removed in its entirety.
If you wish to access the latest MPMP version, you can use
`data("mpmp.28Aug2024")` for your analysis which was scraped by us. If
you wish to use this geneset for MPMP pathway enrichment analysis using
*[pathfindR](https://CRAN.R-project.org/package=pathfindR)*, you can do
so by using `data("pathfindrMPMP")`.

> Note: We urge the users to cite the original articles of the related
> datasets alongside plasmoRUtils.

### Fetching genome metadata and strain names

In above examples, we saw the importance of passing exact name to the
`org` argument for
[`toGeneid()`](https://rohit-satyam.github.io/plasmoRUtils/reference/toGeneid.md)
to function properly. A helper function is provided to achieve this
called
[`listVeupathdb()`](https://rohit-satyam.github.io/plasmoRUtils/reference/listVeupathdb.md).
By default, 11 columns are returned including organism name, the
respective database present in “VEuPathDB Project” column and some more
additional information. However, you can limit the search to columns of
your interests as shown below. As stated by VEuPathDB:

> The best use of the VEuPathDB portal is to get a table with all
> organisms in our sites, and for each organism: the component site, and
> the urls to access their fasta and gff files.

``` r

listVeupathdb() %>% head()
#> # A tibble: 6 × 11
#>   Organism                 Species Genome Fasta Downloa…¹ CDS Fasta Download L…²
#>   <chr>                    <chr>   <chr>                  <chr>                 
#> 1 Acanthamoeba astronyxis… Acanth… http://AmoebaDB.org/a… http://AmoebaDB.org/a…
#> 2 Acanthamoeba castellani… Acanth… http://AmoebaDB.org/a… http://AmoebaDB.org/a…
#> 3 Acanthamoeba castellani… Acanth… http://AmoebaDB.org/a… http://AmoebaDB.org/a…
#> 4 Acanthamoeba castellani… Acanth… http://AmoebaDB.org/a… http://AmoebaDB.org/a…
#> 5 Acanthamoeba castellani… Acanth… http://AmoebaDB.org/a… http://AmoebaDB.org/a…
#> 6 Acanthamoeba culbertson… Acanth… http://AmoebaDB.org/a… http://AmoebaDB.org/a…
#> # ℹ abbreviated names: ¹​`Genome Fasta Download Link`,
#> #   ²​`CDS Fasta Download Link`
#> # ℹ 7 more variables: `Transcript Fasta Download Link` <chr>,
#> #   `Protein Fasta Download Link` <chr>, project_id <chr>, Genes <chr>,
#> #   `GFF Download Link` <chr>, `Genome Source` <chr>,
#> #   `Structural Annotation Source` <chr>
listVeupathdb(customFields=c("species", "project_id")) %>% head()
#> # A tibble: 6 × 3
#>   Organism                                Species                  project_id
#>   <chr>                                   <chr>                    <chr>     
#> 1 Acanthamoeba astronyxis Unknown         Acanthamoeba astronyxis  AmoebaDB  
#> 2 Acanthamoeba castellanii C3             Acanthamoeba castellanii AmoebaDB  
#> 3 Acanthamoeba castellanii Ma             Acanthamoeba castellanii AmoebaDB  
#> 4 Acanthamoeba castellanii str. Neff      Acanthamoeba castellanii AmoebaDB  
#> 5 Acanthamoeba castellanii str. Neff 2021 Acanthamoeba castellanii AmoebaDB  
#> 6 Acanthamoeba culbertsoni A1             Acanthamoeba culbertsoni AmoebaDB
```

### Mapping PDB IDs to Gene IDs

There is currently no facility in VEuPathDB to convert the PDB IDs to
respective gene IDs. If PDB ID corresponds to a multimer complex and if
you have multiple such PDB ids, it becomes arduous to map them to Gene
IDs manually. To provide solution to this issue, we can first convert
the PDB chains to Uniprot IDs using our own
[`pdb2uniprot()`](https://rohit-satyam.github.io/plasmoRUtils/reference/pdb2uniprot.md)
function and then can use
[`toGeneid()`](https://rohit-satyam.github.io/plasmoRUtils/reference/toGeneid.md)
function to obtain gene IDs.

``` r

pdbids <- c("7D2W", "4U5A", "9YGJ")
df <- lapply(pdbids, pdb2uniprot) %>% plyr::ldply()

## Works nicely when all the Uniprot IDs are from same species and strain.
geneids <- toGeneid(inputid = unique(df$attribute),from = "uniprot",to = "ensembl")[,c(1,2)]

## Combining the geneIDs with df
S4Vectors::merge(df,geneids,all=TRUE, by.x="attribute", by.y="Input ID") %>% head()
#>    attribute entity_id chain_id struct_asym_id unp_start unp_end
#> 1 A0A143ZZR8         1        A              A        27     206
#> 2 A0A143ZZR8         1        B              B        27     206
#> 3 A0A2I0BXF6         3        C              C         1    1156
#> 4 A0A2I0C0S7         2        B              B         1    1201
#> 5 A0A2I0C357         1        A              A         1    2189
#> 6     Q75UY1         1        A              A        42     241
#>   start.author_residue_number start.author_insertion_code start.residue_number
#> 1                          NA                                                3
#> 2                          NA                                                3
#> 3                          NA                                                1
#> 4                          NA                                                1
#> 5                          NA                                                1
#> 6                          NA                                                2
#>   end.author_residue_number end.author_insertion_code end.residue_number
#> 1                       181                                          182
#> 2                       181                                          182
#> 3                      1156                                         1156
#> 4                        NA                                         1201
#> 5                        NA                                         2189
#> 6                        NA                                          201
#>   identity coverage query        Output ID
#> 1        1    0.994  7D2W    PF3D7_1372300
#> 2        1    0.994  7D2W    PF3D7_1372300
#> 3        1    1.000  9YGJ PfNF54_080021500
#> 4        1    1.000  9YGJ PfNF54_110020300
#> 5        1    1.000  9YGJ PfNF54_140055700
#> 6        1    0.901  4U5A   PBANKA_1355600
```

## Session

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
#>  [1] tidyselect_1.2.1    viridisLite_0.4.3   farver_2.1.2       
#>  [4] S7_0.2.2            fastmap_1.2.0       mgsub_2.0.0        
#>  [7] promises_1.5.0      digest_0.6.39       timechange_0.4.0   
#> [10] lifecycle_1.0.5     processx_3.9.0      magrittr_2.0.5     
#> [13] compiler_4.4.1      sass_0.4.10         tools_4.4.1        
#> [16] utf8_1.2.6          yaml_2.3.12         gt_1.4.0           
#> [19] data.table_1.18.6.1 knitr_1.52          htmlwidgets_1.6.4  
#> [22] bit_4.6.0           easyPubMed_3.1.6    curl_8.1.0         
#> [25] plyr_1.8.9          xml2_1.6.0          RColorBrewer_1.1-3 
#> [28] websocket_1.4.4     withr_3.0.3         purrr_1.2.2        
#> [31] BiocGenerics_0.50.0 desc_1.4.3          grid_4.4.1         
#> [34] stats4_4.4.1        ggplot2_4.0.3       scales_1.4.0       
#> [37] dichromat_2.0-1     cli_3.6.6           rmarkdown_2.32     
#> [40] crayon_1.5.3        ragg_1.5.2          generics_0.1.4     
#> [43] otel_0.2.0          rstudioapi_0.19.0   httr_1.4.9         
#> [46] reshape2_1.4.5      tzdb_0.5.0          cachem_1.1.0       
#> [49] chromote_0.5.1      stringr_1.6.0       rvest_1.0.5        
#> [52] parallel_4.4.1      BiocManager_1.30.27 vctrs_0.7.3        
#> [55] jsonlite_2.0.0      bookdown_0.48       hms_1.1.4          
#> [58] S4Vectors_0.42.1    bit64_4.8.6         systemfonts_1.3.2  
#> [61] plotly_4.12.1       tidyr_1.3.2         jquerylib_0.1.4    
#> [64] bio3d_2.4-5         glue_1.8.1          pkgdown_2.2.1      
#> [67] codetools_0.2-20    lubridate_1.9.5     stringi_1.8.9      
#> [70] gtable_0.3.6        later_1.4.8         tibble_3.3.1       
#> [73] pillar_1.11.1       htmltools_0.5.9     R6_2.6.1           
#> [76] httr2_1.3.0         textshaping_1.0.5   vroom_1.7.1        
#> [79] evaluate_1.0.5      snakecase_0.11.1    bslib_0.12.0       
#> [82] Rcpp_1.1.2          xfun_0.62           fs_2.1.0           
#> [85] pkgconfig_2.0.3
```

## References

Amos, Beatrice, Cristina Aurrecoechea, Matthieu Barba, et al. 2022.
“VEuPathDB: The Eukaryotic Pathogen, Vector and Host Bioinformatics
Resource Center.” *Nucleic Acids Research* 50 (D1): D898–911.
