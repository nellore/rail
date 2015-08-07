library('derfinder')

## Input options used
maindir
cutoff
readLen
chr

message(Sys.time())
timeinfo <- NULL
timeinfo <- c(timeinfo, list(Sys.time()));

## Load data
load(file.path(maindir, 'CoverageInfo', paste0(chr, 'CovInfo.Rdata')))
fullCov <- list(get(paste0(chr, 'CovInfo')))
names(fullCov) <- chr
timeinfo <- c(timeinfo, list(Sys.time()))
proc.time()
message(Sys.time())

## Identify number of total mapped reads
if(maindir == '/dcs01/ajaffe/Brain/derRuns/railDER/railGEU') {
    bws <- rawFiles(datadir = '/dcs01/ajaffe/Brain/derRuns/railDER/bigwig', samplepatt = 'bw', fileterm = NULL)
} else if (maindir == '/dcs01/ajaffe/Brain/derRuns/railDER/resub') {
    bws <- rawFiles(datadir = '/dcl01/leek/data/geuvadis_rail_v0.1.9/coverage_bigwigs', samplepatt = 'bw', fileterm = NULL)
    bws <- bws[!grepl('mean|median|unique', bws)]
}


names(bws) <- gsub('.bw', '', names(bws))

if(maindir == '/dcs01/ajaffe/Brain/derRuns/railDER/railGEU') {
    counts <- read.table('/dcs01/ajaffe/Brain/derRuns/railDER/all_of_geuvadis_read_counts_v4.4.2015', header = TRUE, sep = '\t')
    mapped <- counts$mapped.read.count[match(names(bws), counts$sample.name)]
    
} else if (maindir == '/dcs01/ajaffe/Brain/derRuns/railDER/resub') {
    counts <- read.table('/dcl01/leek/data/geuvadis_rail_v0.1.9/cross_sample_results/counts.tsv.gz', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    counts$totalMapped <- as.integer(sapply(strsplit(counts$total.mapped.reads, ','), '[[', 1))
    mapped <- counts$totalMapped[match(names(bws), counts$X)]
}

## Match the names
names(mapped) <- names(bws)

## Check distribution
summary(mapped)
summary(mapped) / 1e6

## run regionMatrix
regionMat <- regionMatrix(fullCov, maxClusterGap = 3000L, L = readLen, cutoff = cutoff, returnBP = FALSE, totalMapped = mapped)
timeinfo <- c(timeinfo, list(Sys.time()))

## Save results
save(regionMat, file=paste0('regionMat-cut', cutoff, '-', chr, '.Rdata'))
timeinfo <- c(timeinfo, list(Sys.time()))

## Save time information
save(timeinfo, file=paste0('timeinfo-', chr, '.Rdata'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
