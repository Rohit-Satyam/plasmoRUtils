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

``` r

# Load package and some other useful packages by using
suppressPackageStartupMessages(
  suppressWarnings({
    library(plasmoRUtils)
    library(dplyr)
    library(plyr)}))
```

## Gene ID Conversion

A common challenge in bioinformatics involves mapping between different
identifier systems across databases. For apicomplexan research, this
might include converting between UniProt IDs, legacy gene identifiers,
and current Ensembl gene IDs. The
[`toGeneid()`](https://rohit-satyam.github.io/plasmoRUtils/reference/toGeneid.md)
function addresses this need by retrieving up-to-date annotations from
VEuPathDB databases, supporting bidirectional conversion between various
ID types through flexible parameter specification.

### Retrieving Gene Annotations and Alternative IDs

The
[`toGeneid()`](https://rohit-satyam.github.io/plasmoRUtils/reference/toGeneid.md)
function enables annotation retrieval when provided with Ensembl gene
IDs. By default, it returns essential information including gene names,
symbols, outdated gene IDs, and protein Uniprot IDs. The function’s
versatility extends to supporting custom field requests through its
`customFields` parameter, with available options documented in the
[`getTable()`](https://rohit-satyam.github.io/plasmoRUtils/reference/getTable.md)
help section.

``` r

## Get annotations for list of geneIDs for PF3D7
toGeneid(c("PF3D7_0420300", "PF3D7_0621000"), from="ensembl")
#> # A tibble: 2 × 10
#>   `Gene ID`     `Product Description`        `Gene Strand` `Gene Name or Symbol`
#>   <chr>         <chr>                        <chr>         <chr>                
#> 1 PF3D7_0420300 AP2 domain transcription fa… forward       ApiAP2               
#> 2 PF3D7_0621000 RNA polymerase subunit sigm… forward       ApSigma              
#> # ℹ 6 more variables: `Previous ID(s)` <chr>, `Entrez Gene ID` <chr>,
#> #   `UniProt ID(s)` <chr>, `Protein Length` <chr>, `# TM Domains` <chr>,
#> #   `SignalP Peptide` <chr>
## Get annotations for list of geneIDs for organisms other than PF3D7
toGeneid(inputid = c("TGME49_304740","TGME49_208030"),from="ensembl",org="Toxoplasma gondii ME49", db="toxodb")
#> # A tibble: 2 × 10
#>   `Gene ID`     `Product Description`        `Gene Strand` `Gene Name or Symbol`
#>   <chr>         <chr>                        <chr>         <chr>                
#> 1 TGME49_208030 microneme protein MIC4       forward       MIC4                 
#> 2 TGME49_304740 rhoptry kinase family prote… reverse       ROP35                
#> # ℹ 6 more variables: `Previous ID(s)` <chr>, `Entrez Gene ID` <chr>,
#> #   `UniProt ID(s)` <chr>, `Protein Length` <chr>, `# TM Domains` <chr>,
#> #   `SignalP Peptide` <chr>

## Convert uniprot IDs back to gene IDs. It will also provide Product description and Gene Symbol

toGeneid(inputid = c("Q8I1N6","C6KT48"),from="uniprot",to="ensembl" )
#> # A tibble: 2 × 2
#>   `Gene ID`     `UniProt ID(s)`
#>   <chr>         <chr>          
#> 1 PF3D7_0420300 Q8I1N6         
#> 2 PF3D7_0621000 C6KT48


## Using customFields to get only columns of interest
toGeneid(inputid = c("TGME49_304740","TGME49_208030"),
         from="ensembl",org="Toxoplasma gondii ME49",
         db="toxodb",
         customFields=c("primary_key","predicted_go_component","annotated_go_function"))
#> # A tibble: 2 × 3
#>   `Gene ID`     `Computed GO Components` `Curated GO Functions`
#>   <chr>         <chr>                    <chr>                 
#> 1 TGME49_208030 extracellular region     N/A                   
#> 2 TGME49_304740 N/A                      N/A
```

> **Note:** Successful ID conversion requires precise organism
> nomenclature matching VEuPathDB’s conventions. For example,
> *`Toxoplasma gondii ME49`* must include proper spacing and special
> characters. Invalid query would be *`Toxoplasma gondiiME49`* or short
> forms `TgME49`.

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

For instance, the two protein IDs `Q8I0P6` and `A0A2P1JI60` correspond
to `elongation factor 1-alpha`. Although the genes `PF3D7_1357100` and
`PF3D7_1357000` encode the same protein (thereby same uniprot ID
`Q8I0P6`), they are located on opposite strands: `PF3D7_1357100` is on
the forward strand, whereas `PF3D7_1357000` is on the reverse strand.
They also differ in transcript length. In this example if the task is
Uniprot to Gene ID mapping, users can group the results returned by
[`toGeneid()`](https://rohit-satyam.github.io/plasmoRUtils/reference/toGeneid.md)
function by Uniprot IDs and collapse rows thereby concatinating two gene
IDs separated by semicolon or comma. An alternative is to randomly
choose one lexicographically. Because there is no universally
appropriate rule for resolving this type of ambiguity, we avoid making
this choice for the end user.

``` r

toGeneid(c("Q8I0P6" ,"A0A2P1JI60"), from = "uniprot",to = "ensembl")
#> # A tibble: 3 × 2
#>   `Gene ID`     `UniProt ID(s)`
#>   <chr>         <chr>          
#> 1 PF3D7_1357000 Q8I0P6         
#> 2 PF3D7_1357100 A0A2P1JI60     
#> 3 PF3D7_1357100 Q8I0P6
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

getTable(org="Plasmodium falciparum 3D7", db="plasmodb") %>% head()
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
#> #   `SignalP Peptide` <chr>

## User can also provide custom fields. For example we wish to download the P. falciparum 3D7 Proteome and phosphoproteome data during intraerythrocytic development (Quantitative) (Pease et al.)

getTable(org="Plasmodium falciparum 3D7", db="plasmodb", customFields = c("primary_key","pan_6365","pan_6366","pan_6367")) %>% head()
#> # A tibble: 6 × 4
#>   `Gene ID` Ring Ave (Global pro…¹ Troph Ave (Global pr…² Schizont Ave (Global…³
#>   <chr>     <chr>                  <chr>                  <chr>                 
#> 1 PF3D7_01… N/A                    N/A                    N/A                   
#> 2 PF3D7_01… N/A                    N/A                    N/A                   
#> 3 PF3D7_01… N/A                    N/A                    N/A                   
#> 4 PF3D7_01… N/A                    N/A                    N/A                   
#> 5 PF3D7_01… N/A                    N/A                    N/A                   
#> 6 PF3D7_01… N/A                    N/A                    N/A                   
#> # ℹ abbreviated names: ¹​`Ring Ave (Global proteome and phosphoproteome)`,
#> #   ²​`Troph Ave (Global proteome and phosphoproteome)`,
#> #   ³​`Schizont Ave (Global proteome and phosphoproteome)`
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

The constituent databases of VEuPathDB also provide some preconfigured
tables which can not be fetched via
[`getTable()`](https://rohit-satyam.github.io/plasmoRUtils/reference/getTable.md)
function. To enable users to fetch such tables, we wrote another
function called
[`getPreconfiguredTable()`](https://rohit-satyam.github.io/plasmoRUtils/reference/getPreconfiguredTable.md)
the usage of which has been shown below.

``` r

## Fetch pathway table for all the genes from MPMP database

getPreconfiguredTable(org = "Plasmodium falciparum 3D7",db = "plasmodb",customField = "MetabolicPathwaysMPMP") %>% head()
#> # A tibble: 6 × 4
#>   `Gene ID`     pathway_id Pathway                                      Activity
#>   <chr>         <chr>      <chr>                                        <chr>   
#> 1 PF3D7_0100100 lys_met    Peptides with confirmed methylated lysine r… erythro…
#> 2 PF3D7_0100100 virulence  Candidate genes related to virulence         erythro…
#> 3 PF3D7_0100100 Par_RBC    Protein-Protein Interactions between Human … erythro…
#> 4 PF3D7_0100100 PfEMP1     PfEMP1 domain architectures                  erythro…
#> 5 PF3D7_0100100 gene_lumef Gene expression affected by lumefantrine     erythro…
#> 6 PF3D7_0100100 PQS        P. falciparum genes harboring G-quadruplexes erythro…
```

Please note that the MPMP pathway version provided by PlasmoDB is
[outdated](https://plasmodb.org/plasmo/app/record/dataset/DS_1e177b728b)
(`03-2019`). Some pathways have been revised or removed in its entirety.
If you wish to access the latest MPMP version, you can use
`data("mpmp.28Aug2024")` for you analysis which was scraped by us. If
you wish to use this geneset for MPMP pathway enrichment analysis using
*[pathfindR](https://CRAN.R-project.org/package=pathfindR)*, you can do
so by using `data("pathfindrMPMP")`.

Similarly, predictions like `TMHMM` and `SignalP` and `InterPro` have
not been updated given the recent funding crunch and should be therefore
used with caution. We will discuss more about it in a separate tutorial.

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
#> 1 Edhazardia aedis USNM 4… Edhaza… http://MicrosporidiaD… http://MicrosporidiaD…
#> 2 Kluyveromyces marxianus… Kluyve… http://FungiDB.org/co… http://FungiDB.org/co…
#> 3 Aspergillus luchuensis … Asperg… http://FungiDB.org/co… http://FungiDB.org/co…
#> 4 Epichloe glyceriae E277  Epichl… http://FungiDB.org/co… http://FungiDB.org/co…
#> 5 Aspergillus versicolor … Asperg… http://FungiDB.org/co… http://FungiDB.org/co…
#> 6 Aspergillus sydowii CBS… Asperg… http://FungiDB.org/co… http://FungiDB.org/co…
#> # ℹ abbreviated names: ¹​`Genome Fasta Download Link`,
#> #   ²​`CDS Fasta Download Link`
#> # ℹ 7 more variables: `Transcript Fasta Download Link` <chr>,
#> #   `Protein Fasta Download Link` <chr>, `VEuPathDB Project` <chr>,
#> #   Genes <chr>, `GFF Download Link` <chr>, `Genome Source` <chr>,
#> #   `Structural Annotation Source` <chr>
listVeupathdb(customFields=c("species", "project_id")) %>% head()
#> # A tibble: 6 × 2
#>   Species                  `VEuPathDB Project`
#>   <chr>                    <chr>              
#> 1 Entamoeba nuttalli       AmoebaDB           
#> 2 Acanthamoeba castellanii AmoebaDB           
#> 3 Mastigamoeba balamuthi   AmoebaDB           
#> 4 Acanthamoeba sp.         AmoebaDB           
#> 5 Acanthamoeba sp.         AmoebaDB           
#> 6 Acanthamoeba sp.         AmoebaDB
```

Since this function also provide URLs of FASTA and GFF files, you can
use it to find the URLs of the files you are interested in and import
them in R directly without leaving the console.

``` r

listVeupathdb() %>% 
  subset(.,Organism =="Edhazardia aedis USNM 41457") %>% 
  select(`GFF Download Link`) %>% as.character() %>% 
  rtracklayer::import.gff3() %>% head()
#> GRanges object with 6 ranges and 12 metadata columns:
#>            seqnames    ranges strand |    source                type     score
#>               <Rle> <IRanges>  <Rle> |  <factor>            <factor> <numeric>
#>   [1] AFBI030000... 5287-6076      + | VEuPathDB protein_coding_gene        NA
#>   [2] AFBI030000... 5287-6076      + | VEuPathDB mRNA                       NA
#>   [3] AFBI030000... 5287-6076      + | VEuPathDB exon                       NA
#>   [4] AFBI030000... 5447-6043      + | VEuPathDB CDS                        NA
#>   [5] AFBI030000... 5287-5446      + | VEuPathDB five_prime_UTR             NA
#>   [6] AFBI030000... 6044-6076      + | VEuPathDB three_prime_UTR            NA
#>           phase            ID   description   ebi_biotype          Parent
#>       <integer>   <character>   <character>   <character> <CharacterList>
#>   [1]      <NA>    EDEG_00001 hypothetic... protein_co...                
#>   [2]      <NA> EDEG_00001... hypothetic...          <NA>      EDEG_00001
#>   [3]      <NA> exon_EDEG_...          <NA>          <NA>   EDEG_00001...
#>   [4]         0 EDEG_00001...          <NA>          <NA>   EDEG_00001...
#>   [5]      <NA> utr_EDEG_0...          <NA>          <NA>   EDEG_00001...
#>   [6]      <NA> utr_EDEG_0...          <NA>          <NA>   EDEG_00001...
#>       gene_ebi_biotype     gene_id protein_source_id            Note
#>            <character> <character>       <character> <CharacterList>
#>   [1]             <NA>        <NA>              <NA>                
#>   [2]    protein_co...        <NA>              <NA>                
#>   [3]             <NA>  EDEG_00001              <NA>                
#>   [4]             <NA>  EDEG_00001     EDEG_00001...                
#>   [5]             <NA>        <NA>              <NA>                
#>   [6]             <NA>        <NA>              <NA>                
#>   -------
#>   seqinfo: 342 sequences from an unspecified genome; no seqlengths
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

pdbids <- c("7D2W","4U5A","6E10")
df <- lapply(pdbids, pdb2uniprot) %>% plyr::ldply()

geneids <- toGeneid(inputid = unique(df$attribute),from = "uniprot",to = "ensembl")

## Combining the geneIDs with df
S4Vectors::merge(df,geneids,all=TRUE, by.x="attribute", by.y="UniProt ID(s)") %>% head()
#>    attribute entity_id chain_id struct_asym_id unp_start unp_end
#> 1 A0A143ZZR8         1        A              A        27     206
#> 2 A0A143ZZR8         1        B              B        27     206
#> 3     Q75UY1         1        C              C        42     241
#> 4     Q75UY1         1        D              D        42     241
#> 5     Q75UY1         1        A              A        42     241
#> 6     Q75UY1         1        B              B        42     241
#>   start.author_residue_number start.author_insertion_code start.residue_number
#> 1                          NA                                                3
#> 2                          NA                                                3
#> 3                          NA                                                2
#> 4                          NA                                                2
#> 5                          NA                                                2
#> 6                          NA                                                2
#>   end.author_residue_number end.author_insertion_code end.residue_number
#> 1                       181                                          182
#> 2                       181                                          182
#> 3                        NA                                          201
#> 4                        NA                                          201
#> 5                        NA                                          201
#> 6                        NA                                          201
#>   identity coverage query       Gene ID
#> 1        1    0.994  7D2W PF3D7_1372300
#> 2        1    0.994  7D2W PF3D7_1372300
#> 3        1    0.901  4U5A          <NA>
#> 4        1    0.901  4U5A          <NA>
#> 5        1    0.901  4U5A          <NA>
#> 6        1    0.901  4U5A          <NA>
```

## Session

``` r

sessionInfo()
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
#> [1] plyr_1.8.9         dplyr_1.2.1        plasmoRUtils_1.1.0 rlang_1.3.0       
#> [5] readr_2.2.0        janitor_2.2.1      BiocStyle_2.32.1  
#> 
#> loaded via a namespace (and not attached):
#>   [1] segmented_2.2-1             fs_2.1.0                   
#>   [3] ProtGenerics_1.36.0         matrixStats_1.5.0          
#>   [5] bitops_1.0-9                lubridate_1.9.5            
#>   [7] pRoloc_1.44.1               httr_1.4.8                 
#>   [9] RColorBrewer_1.1-3          doParallel_1.0.17          
#>  [11] ggsci_5.1.0                 tools_4.4.1                
#>  [13] MSnbase_2.30.1              backports_1.5.1            
#>  [15] utf8_1.2.6                  R6_2.6.1                   
#>  [17] lazyeval_0.2.3              withr_3.0.3                
#>  [19] prettyunits_1.2.0           gridExtra_2.3.1            
#>  [21] preprocessCore_1.66.0       cli_3.6.6                  
#>  [23] Biobase_2.64.0              textshaping_1.0.5          
#>  [25] gt_1.3.0                    sass_0.4.10                
#>  [27] topGO_2.56.0                mvtnorm_1.4-1              
#>  [29] S7_0.2.2                    randomForest_4.7-1.2       
#>  [31] proxy_0.4-29                pkgdown_2.2.0              
#>  [33] Rsamtools_2.20.0            systemfonts_1.3.2          
#>  [35] txdbmaker_1.0.1             AnnotationForge_1.46.0     
#>  [37] dichromat_2.0-0.1           parallelly_1.48.0          
#>  [39] limma_3.60.6                rstudioapi_0.19.0          
#>  [41] impute_1.78.0               RSQLite_3.53.3             
#>  [43] FNN_1.1.4.1                 generics_0.1.4             
#>  [45] BiocIO_1.14.0               vroom_1.7.1                
#>  [47] gtools_3.9.5                car_3.1-5                  
#>  [49] dendextend_1.19.1           GO.db_3.19.1               
#>  [51] Matrix_1.7-5                MALDIquant_1.22.3          
#>  [53] drawProteins_1.24.0         S4Vectors_0.42.1           
#>  [55] abind_1.4-8                 lifecycle_1.0.5            
#>  [57] yaml_2.3.12                 snakecase_0.11.1           
#>  [59] carData_3.0-6               SummarizedExperiment_1.34.0
#>  [61] recipes_1.3.3               SparseArray_1.4.8          
#>  [63] BiocFileCache_2.12.0        grid_4.4.1                 
#>  [65] blob_1.3.0                  promises_1.5.0             
#>  [67] crayon_1.5.3                PSMatch_1.8.0              
#>  [69] lattice_0.22-9              beachmat_2.20.0            
#>  [71] annotate_1.82.0             GenomicFeatures_1.56.0     
#>  [73] chromote_0.5.1              mzR_2.38.0                 
#>  [75] KEGGREST_1.44.1             pillar_1.11.1              
#>  [77] knitr_1.51                  GenomicRanges_1.56.2       
#>  [79] rjson_0.2.23                lpSolve_5.6.23             
#>  [81] future.apply_1.20.2         codetools_0.2-20           
#>  [83] mgsub_2.0.0                 glue_1.8.1                 
#>  [85] pcaMethods_1.96.0           data.table_1.18.4          
#>  [87] MultiAssayExperiment_1.30.3 vctrs_0.7.3                
#>  [89] png_0.1-9                   gtable_0.3.6               
#>  [91] kernlab_0.9-33              cachem_1.1.0               
#>  [93] gower_1.0.2                 xfun_0.59                  
#>  [95] prodlim_2026.03.11          S4Arrays_1.4.1             
#>  [97] polyglotr_1.7.4             coda_0.19-4.1              
#>  [99] survival_3.8-6              ncdf4_1.24                 
#> [101] timeDate_4052.112           SingleCellExperiment_1.26.0
#> [103] iterators_1.0.14            hardhat_1.4.3              
#> [105] lava_1.9.2                  statmod_1.5.2              
#> [107] MLInterfaces_1.84.0         ipred_0.9-15               
#> [109] nlme_3.1-169                bit64_4.8.2                
#> [111] progress_1.2.3              filelock_1.0.3             
#> [113] LaplacesDemon_16.1.8        GenomeInfoDb_1.40.1        
#> [115] bslib_0.11.0                affyio_1.74.0              
#> [117] irlba_2.3.7                 rpart_4.1.27               
#> [119] otel_0.2.0                  colorspace_2.1-2           
#> [121] BiocGenerics_0.50.0         DBI_1.3.0                  
#> [123] nnet_7.3-20                 tidyselect_1.2.1           
#> [125] processx_3.9.0              bit_4.6.0                  
#> [127] compiler_4.4.1              curl_7.1.0                 
#> [129] rvest_1.0.5                 httr2_1.2.3                
#> [131] graph_1.82.0                SparseM_1.84-2             
#> [133] xml2_1.6.0                  desc_1.4.3                 
#> [135] DelayedArray_0.30.1         plotly_4.12.0              
#> [137] bookdown_0.47               rtracklayer_1.64.0         
#> [139] scales_1.4.0                hexbin_1.28.5              
#> [141] affy_1.82.0                 rappdirs_0.3.4             
#> [143] stringr_1.6.0               digest_0.6.39              
#> [145] mixtools_2.0.0.1            rmarkdown_2.31             
#> [147] XVector_0.44.0              htmltools_0.5.9            
#> [149] pkgconfig_2.0.3             SingleR_2.6.0              
#> [151] sparseMatrixStats_1.16.0    MatrixGenerics_1.16.0      
#> [153] dbplyr_2.6.0                fastmap_1.2.0              
#> [155] htmlwidgets_1.6.4           UCSC.utils_1.0.0           
#> [157] DelayedMatrixStats_1.26.0   farver_2.1.2               
#> [159] jquerylib_0.1.4             jsonlite_2.0.0             
#> [161] mclust_6.1.3                BiocParallel_1.38.0        
#> [163] mzID_1.42.0                 ModelMetrics_1.2.2.2       
#> [165] BiocSingular_1.20.0         RCurl_1.98-1.19            
#> [167] magrittr_2.0.5              scuttle_1.14.0             
#> [169] Formula_1.2-5               GenomeInfoDbData_1.2.12    
#> [171] Rcpp_1.1.2                  viridis_0.6.5              
#> [173] MsCoreUtils_1.16.1          vsn_3.72.0                 
#> [175] pROC_1.19.0.1               stringi_1.8.7              
#> [177] zlibbioc_1.50.0             MASS_7.3-65                
#> [179] listenv_1.0.0               parallel_4.4.1             
#> [181] Biostrings_2.72.1           splines_4.4.1              
#> [183] hms_1.1.4                   igraph_2.3.3               
#> [185] ggpubr_1.0.0                QFeatures_1.14.2           
#> [187] ggsignif_0.6.4              reshape2_1.4.5             
#> [189] biomaRt_2.60.1              stats4_4.4.1               
#> [191] ScaledMatrix_1.12.0         XML_3.99-0.23              
#> [193] evaluate_1.0.5              BiocManager_1.30.27        
#> [195] tzdb_0.5.0                  foreach_1.5.2              
#> [197] tidyr_1.3.2                 purrr_1.2.2                
#> [199] future_1.70.0               clue_0.3-68                
#> [201] bio3d_2.4-5                 ggplot2_4.0.3              
#> [203] rsvd_1.0.5                  xtable_1.8-8               
#> [205] broom_1.0.13                restfulr_0.0.17            
#> [207] AnnotationFilter_1.28.0     easyPubMed_3.1.6           
#> [209] e1071_1.7-17                rstatix_1.0.0              
#> [211] later_1.4.8                 class_7.3-23               
#> [213] viridisLite_0.4.3           ragg_1.5.2                 
#> [215] tibble_3.3.1                websocket_1.4.4            
#> [217] memoise_2.0.1               AnnotationDbi_1.66.0       
#> [219] GenomicAlignments_1.40.0    IRanges_2.38.1             
#> [221] cluster_2.1.8.2             globals_0.19.1             
#> [223] timechange_0.4.0            caret_7.0-1                
#> [225] sampling_2.11
```

## References

Amos, Beatrice, Cristina Aurrecoechea, Matthieu Barba, et al. 2022.
“VEuPathDB: The Eukaryotic Pathogen, Vector and Host Bioinformatics
Resource Center.” *Nucleic Acids Research* 50 (D1): D898–911.
