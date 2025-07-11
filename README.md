
<!-- README.md is generated from README.Rmd. Please edit that file -->

# plasmoRUtils <img  src="man/figures/logo.png" align="right" geight="139"/>

`plasmoRUtils` enables users to connect to several *Plasmodium* and
other apicomplexan databases via R interface and provides simple
functions to carry out other bioinformatics tasks which are non-trival
for parasite informatic analysis. For further details, we recommend you
read our preprint.

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![codecov](https://codecov.io/gh/Rohit-Satyam/plasmoRUtils/branch/master/graph/badge.svg)](https://codecov.io/gh/Rohit-Satyam/plasmoRUtils)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://choosealicense.com/licenses/mit/)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-4-6666ff.svg)](https://cran.r-project.org/)
[![packageversion](https://img.shields.io/badge/Package%20version-1.0.0-turquoise.svg?style=flat-square)](commits/master)
[![Last-changedate](https://img.shields.io/badge/last%20change-2025--06--30-yellowgreen.svg)](/commits/master)
<!-- badges: end -->

## Installation

Before downloading the package, install the following dependencies.

``` r
cranpkgs <- c('BiocManager','randomcoloR', 'janitor', 'readr', 'rlang', 'dplyr', 'ggsci', 'rvest',
'easyPubMed', 'plyr', 'scales', 'ggplot2', 'glue', 'tidyr', 'tibble', 'data.table', 'plotly',
'purrr', 'stringr', 'S4Vectors', 'echarts4r', 'magrittr', 'bio3d', 'httr', 'jsonlite',
'ggpubr', 'gt', 'mgsub', 'reshape2','pathfindR')

install.packages(setdiff(cranpkgs, rownames(installed.packages())), dependencies = TRUE)

biocpkgs <- c("rmarkdown","pRoloc","knitr","BiocStyle","DESeq2","styler","utils","IRanges",
"BiocGenerics","rtracklayer","scuttle","txdbmaker","topGO","drawProteins","GenomicFeatures",
"biomaRt","AnnotationForge","Biostrings","GenomeInfoDb","SingleCellExperiment",
"SingleR","NOISeq","GenomicRanges","BSgenome")

BiocManager::install(setdiff(biocpkgs, rownames(installed.packages())), dependencies = TRUE)
```

You can install the development version of `plasmoRUtils` using:

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
The detailed vignettes are available at `plasmoRUtils` website. To visit, click [here](https://rohit-satyam.github.io/plasmoRUtils/)

The documentation of this package is also available at the following:

1.  [Introduction to
    plasmoRUtils](https://htmlpreview.github.io/?https://github.com/Rohit-Satyam/plasmoRUtils/main/vignettes/Introduction_to_plasmoRUtils.html)
2.  [Accessing component databases of
    VEuPathDB](https://htmlpreview.github.io/?https://github.com/Rohit-Satyam/plasmoRUtils/main/vignettes/Gene_ID_Conversion.html)
3.  [Other useful
    functions](https://htmlpreview.github.io/?https://github.com/Rohit-Satyam/plasmoRUtils/main/vignettes/Miscellaneous_function.html)

## To-do List

1.  Provide a function to access Plasmobase.
2.  Write a wrapper function, easypathFindR, to perform Pathway
    enrichment analysis quickly.
3.  Write a function to make String PPI quickly.

## Contributing

We’re excited to have you contribute to this package! If you’d like to
help out, try to follow the same style and conventions used in the
current functions - where it makes sense, of course. If you have any
ideas or suggestions, don’t hesitate to reach out—opening a GitHub issue
is usually the best way to start the conversation.

> Just a heads up: this project has a Contributor Code of Conduct, so by
> getting involved, you’re agreeing to play by those rules. Thanks for
> helping make this project better!
