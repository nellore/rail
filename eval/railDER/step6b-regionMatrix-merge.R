library('derfinder')

## Input options used
cutoff

## Load data
chrs <- paste0('chr', c(1:22, 'X', 'Y'))
regionMat <- lapply(chrs, function(chr) {
    load(paste0('regionMat-cut', cutoff, '-', chr, '.Rdata'))
    res <- regionMat
    return(res)
})

## Merge
regionMat <- do.call(c, regionMat)

## Save
save(regionMat, file = paste0('regionMat-cut', cutoff, '.Rdata'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
