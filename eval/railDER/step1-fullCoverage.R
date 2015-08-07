## Load the data without a filter, save it, then filter it for derfinder processing steps

## Load libraries
library('getopt')

## Available at http://www.bioconductor.org/packages/release/bioc/html/derfinder.html
library('derfinder')
library('BiocParallel')
library('devtools')
library('rtracklayer')

## Specify parameters
spec <- matrix(c(
	'datadir', 'd', 1, 'character', 'Data directory, matched with rawFiles(datadir)',
	'pattern', 'p', 1, 'character', 'Sample pattern',
	'cutoff', 'c', 1, 'numeric', 'Filtering cutoff used',
	'mcores', 'm', 1, 'integer', 'Number of cores',
    'fileStyle', 'f', 2, 'character', 'FileStyle used for naming the chromosomes',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## Default value for fileStyle
if (is.null(opt$fileStyle)) opt$fileStyle <- 'UCSC'

## Identify the data files
bws <- rawFiles(datadir = opt$datadir, samplepatt = opt$pattern, fileterm = NULL)
bws <- bws[!grepl('mean|median|unique', bws)]
names(bws) <- gsub('.bw', '', names(bws))

bws <- BigWigFileList(bws)

## Load the coverage information without filtering
chrnums <- c(1:22, 'X', 'Y')

fullCov <- fullCoverage(files = bws, chrs = chrnums, mc.cores = opt$mcores, fileStyle = opt$fileStyle, outputs = 'auto')

message(paste(Sys.time(), 'Saving the full (unfiltered) coverage data'))
save(fullCov, file='fullCov.Rdata')

## Filter the data and save it by chr
myFilt <- function(chr, rawData, cutoff) {
    library('derfinder')
    message(paste(Sys.time(), 'Filtering chromosome', chr))
    
	## Filter the data
	res <- filterData(data = rawData, cutoff = cutoff, index = NULL)
	
	## Save it in a unified name format
	varname <- paste0(chr, 'CovInfo')
	assign(varname, res)
	output <- paste0(varname, '.Rdata')
	
	## Save the filtered data
	save(list = varname, file = output, compress='gzip')
	
	## Finish
	return(invisible(NULL))
}

#message(paste(Sys.time(), 'Filtering and saving the data with cutoff', opt$cutoff))
#filteredCov <- bpmapply(myFilt, names(fullCov), fullCov, BPPARAM = SnowParam(opt$mcores, outfile = Sys.getenv('SGE_STDERR_PATH')), MoreArgs = list(cutoff = opt$cutoff))

## Done!
proc.time()
options(width = 120)
session_info()
