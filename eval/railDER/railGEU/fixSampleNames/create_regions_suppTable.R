library('GenomicRanges')
library('devtools')
load('regions.Rdata')
names(regions) <- NULL
write.csv(as.data.frame(regions), file = "supplementaryExpressedRegions.csv", 
    quote = FALSE, row.names = FALSE)
options(width = 120)
session_info()

#> setwd('/Users/lcollado/Dropbox/JHSPH/Code/rail/eval/railDER/railGEU/fixSampleNames')
#> library('GenomicRanges')
#Loading required package: BiocGenerics
#Loading required package: parallel
#
#Attaching package: ‘BiocGenerics’
#
#The following objects are masked from ‘package:parallel’:
#
#    clusterApply, clusterApplyLB, clusterCall, clusterEvalQ, clusterExport, clusterMap, parApply, parCapply, parLapply, parLapplyLB,
#    parRapply, parSapply, parSapplyLB
#
#The following object is masked from ‘package:stats’:
#
#    xtabs
#
#The following objects are masked from ‘package:base’:
#
#    anyDuplicated, append, as.data.frame, as.vector, cbind, colnames, do.call, duplicated, eval, evalq, Filter, Find, get, intersect,
#    is.unsorted, lapply, Map, mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce,
#    rep.int, rownames, sapply, setdiff, sort, table, tapply, union, unique, unlist, unsplit
#
#Loading required package: S4Vectors
#Loading required package: stats4
#Loading required package: IRanges
#Loading required package: GenomeInfoDb
#> library('devtools')
#> load('regions.Rdata')
#> names(regions) <- NULL
#> write.csv(as.data.frame(regions), file = "supplementaryExpressedRegions.csv", 
#+     quote = FALSE, row.names = FALSE)
#> options(width = 120)
#> session_info()
#Session info-----------------------------------------------------------------------------------------------------------
# setting  value                                             
# version  R Under development (unstable) (2014-11-01 r66923)
# system   x86_64, darwin10.8.0                              
# ui       AQUA                                              
# language (EN)                                              
# collate  en_US.UTF-8                                       
# tz       America/Detroit                                   
#
#Packages---------------------------------------------------------------------------------------------------------------
# package       * version date       source        
# BiocGenerics  * 0.13.11 2015-04-03 Bioconductor  
# devtools      * 1.6.1   2014-10-07 CRAN (R 3.2.0)
# GenomeInfoDb  * 1.3.16  2015-03-27 Bioconductor  
# GenomicRanges * 1.19.52 2015-04-04 Bioconductor  
# IRanges       * 2.1.43  2015-03-07 Bioconductor  
# rstudioapi      0.3.1   2015-04-07 CRAN (R 3.2.0)
# S4Vectors     * 0.5.22  2015-03-06 Bioconductor  
# XVector         0.7.4   2015-02-08 Bioconductor
