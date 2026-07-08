# plasmoRUtils

`plasmoRUtils` enables users to connect to several *Plasmodium* and
other apicomplexan databases via R interface and provides simple
functions to carry out other bioinformatics tasks which are non-trival
for parasite bioinformatic analysis. For further details, we recommend
you read our
[preprint](https://www.biorxiv.org/content/10.1101/2025.07.30.667718v1).

## Installation

The easiest way to download package is as follows

``` r

pak::pkg_install("Rohit-Satyam/plasmoRUtils", dependencies = TRUE)
```

If the above method fails, try the following steps:

1.  Before downloading the package, install the following dependencies.

``` r

cranpkgs <- c('BiocManager','randomcoloR', 'janitor', 'readr', 'rlang', 'dplyr', 'rvest',
'easyPubMed', 'plyr', 'scales', 'ggplot2', 'tidyr', 'tibble', 'data.table', 'plotly',
'purrr', 'stringr', 'S4Vectors', 'magrittr', 'bio3d', 'httr2', 'jsonlite', 'gt', 'mgsub', 'reshape2','pathfindR')

install.packages(setdiff(cranpkgs, rownames(installed.packages())), dependencies = TRUE)

biocpkgs <- c("rmarkdown","pRoloc","knitr","BiocStyle","DESeq2","styler","utils","IRanges",
"BiocGenerics","rtracklayer","scuttle","txdbmaker","topGO","drawProteins","GenomicFeatures",
"biomaRt","AnnotationForge","Biostrings","GenomeInfoDb","SingleCellExperiment",
"SingleR","NOISeq","GenomicRanges","BSgenome")

BiocManager::install(setdiff(biocpkgs, rownames(installed.packages())), dependencies = TRUE)
```

2.  You can install the development version of `plasmoRUtils` using:

``` r

devtools::install_github("Rohit-Satyam/plasmoRUtils")
remotes::install_github('Rohit-Satyam/plasmoRUtils')
```

## Check installation

Once dependencies are installed, the package can be loaded as follows:

``` r

# Once installed load the library as
library(plasmoRUtils)

## To re-check if all the dependencies that are required by plasmoRUtils are installed
install_dependencies()
```

## Documentation

The documentation of this package is available at the following:

1.  [Introduction to
    plasmoRUtils](https://rohit-satyam.github.io/plasmoRUtils/articles/Introduction_to_plasmoRUtils.html)
2.  [Accessing component databases of
    VEuPathDB](https://rohit-satyam.github.io/plasmoRUtils/articles/Gene_ID_Conversion.html)
3.  [Other useful
    functions](https://rohit-satyam.github.io/plasmoRUtils/articles/Miscellaneous_function.html)
4.  [RNASeq: Importance of
    reanalysis](https://rohit-satyam.github.io/plasmoRUtils/articles/Need_for_reanalysis.html)
5.  [Microarray data
    reanalysis](https://rohit-satyam.github.io/plasmoRUtils/articles/Microarray_reanalysis.html)

## To-do List

1.  Write a wrapper function, easypathFindR, to perform Pathway
    enrichment analysis quickly.
2.  Write a function to make String PPI quickly.

## Contributing

We’re excited to have you contribute to this package! If you’d like to
help out, try to follow the same style and conventions used in the
current functions - where it makes sense, of course. If you have any
ideas or suggestions, don’t hesitate to reach out—opening a GitHub issue
is usually the best way to start the conversation.

> Just a heads up: this project has a Contributor Code of Conduct, so by
> getting involved, you’re agreeing to play by those rules. Thanks for
> helping make this project better!
