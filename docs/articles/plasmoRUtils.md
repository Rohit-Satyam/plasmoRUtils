# plasmoRUtils

``` r

library(plasmoRUtils)
#> Loading required package: pRoloc
#> Loading required package: MSnbase
#> Loading required package: BiocGenerics
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
#>     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     Position, rank, rbind, Reduce, rownames, sapply, setdiff, table,
#>     tapply, union, unique, unsplit, which.max, which.min
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> Loading required package: mzR
#> Loading required package: Rcpp
#> Warning in fun(libname, pkgname): mzR has been built against a different Rcpp version (1.0.12)
#> than is installed on your system (1.0.13). This might lead to errors
#> when loading mzR. If you encounter such issues, please send a report,
#> including the output of sessionInfo() to the Bioc support forum at 
#> https://support.bioconductor.org/. For details see also
#> https://github.com/sneumann/mzR/wiki/mzR-Rcpp-compiler-linker-issue.
#> Loading required package: S4Vectors
#> Loading required package: stats4
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: ProtGenerics
#> 
#> Attaching package: 'ProtGenerics'
#> The following object is masked from 'package:stats':
#> 
#>     smooth
#> 
#> This is MSnbase version 2.30.1 
#>   Visit https://lgatto.github.io/MSnbase/ to get started.
#>  Consider switching to the 'R for Mass Spectrometry'
#>  packages - see https://RforMassSpectrometry.org for details.
#> 
#> Attaching package: 'MSnbase'
#> The following object is masked from 'package:base':
#> 
#>     trimws
#> Loading required package: MLInterfaces
#> Loading required package: annotate
#> Loading required package: AnnotationDbi
#> Loading required package: IRanges
#> 
#> Attaching package: 'IRanges'
#> The following object is masked from 'package:grDevices':
#> 
#>     windows
#> Loading required package: XML
#> 
#> Attaching package: 'annotate'
#> The following object is masked from 'package:mzR':
#> 
#>     nChrom
#> Loading required package: cluster
#> Loading required package: BiocParallel
#> 
#> This is pRoloc version 1.44.1 
#>   Visit https://lgatto.github.io/pRoloc/ to get started.
#> Loading required package: janitor
#> 
#> Attaching package: 'janitor'
#> The following objects are masked from 'package:stats':
#> 
#>     chisq.test, fisher.test
#> Loading required package: readr
#> Loading required package: rlang
#> 
#> Attaching package: 'rlang'
#> The following object is masked from 'package:MSnbase':
#> 
#>     exprs
#> The following object is masked from 'package:Biobase':
#> 
#>     exprs
#> Warning: replacing previous import 'data.table::first' by 'dplyr::first' when
#> loading 'plasmoRUtils'
#> Warning: replacing previous import 'biomaRt::select' by 'dplyr::select' when
#> loading 'plasmoRUtils'
#> Warning: replacing previous import 'data.table::last' by 'dplyr::last' when
#> loading 'plasmoRUtils'
#> Warning: replacing previous import 'data.table::between' by 'dplyr::between'
#> when loading 'plasmoRUtils'
#> Warning: replacing previous import 'data.table::transpose' by
#> 'purrr::transpose' when loading 'plasmoRUtils'
#> 
#> Warning: replacing previous import 'GenomicFeatures::genes' by 'topGO::genes'
#> when loading 'plasmoRUtils'
```
