source("http://bioconductor.org/biocLite.R")
biocLite("rtracklayer")
library(rtracklayer)
r1 <- IRanges::RangesList(chr1 = IRanges::IRanges(c(1, 5), c(3, 6)))
a <- import("http://burn.fm/bigwigs/PRJEB3366_ERP001942_ERS185031_ERX162976_ERR188027-1-1.bw", format="bw", selection=BigWigSelection(r1))
a