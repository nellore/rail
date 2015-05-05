library('GenomicRanges')
library('ballgown')
load('../regionMatrix/regionMat-cut5.Rdata')
# Available from https://www.dropbox.com/s/kp5th9hgkq8ckom/geuvadisbg.rda
load('geuvadisbg.rda')

## Pheno data
pd <- pData(geuvadisbg)

## Extract all regions
regions <- unlist(GRangesList(lapply(regionMat, '[[', 'regions')))

## Extract all coverageMatrices
coverageMatrix <- do.call(rbind, lapply(regionMat, '[[', 'coverageMatrix'))
colnames(coverageMatrix) <- gsub('\\.', '-', colnames(coverageMatrix))
cNames <- colnames(coverageMatrix)

## Load https://github.com/buci/rail/blob/master/eval/E-GEUV-3.sdrf.txt
sdrf <- read.table('E-GEUV-3.sdrf.txt', header = TRUE, sep = '\t')
sdrf <- subset(sdrf, Derived.Array.Data.File != "")
bam <- gsub('.bam', '', sdrf$Derived.Array.Data.File)

## Load matching info from manifest file
## https://github.com/buci/rail/blob/master/eval/GEUVADIS_all_descriptive.manifest
manifest <- read.table('GEUVADIS_all_descriptive.manifest', skip = 6, sep = '\t')

## Match sample names with those from the ballgown object: rail to bg
match_c_man <- match(cNames, manifest$V5)
match_c_sd <- match(manifest$V1[match_c_man], sdrf$Comment.FASTQ_URI.)
match_c_bg <- match(bam[match_c_sd], pd$SampleID)

## Reverse match: bg to rail
match_bg_sd <- match(pd$SampleID, bam)
match_bg_man <- match(sdrf$Comment.FASTQ_URI.[match_bg_sd], manifest$V1)
match_bg_c <- match(manifest$V5[match_bg_man], cNames)

## Save matched info
pMatch <- data.frame('railName' = cNames, 'bgName' = pd$SampleID[match_c_bg])

## Sort cov by bgName
coverageMatrix <- coverageMatrix[, match_bg_c]

save(pMatch, file = 'pMatch.Rdata')
save(match_bg_c, match_c_bg, file = 'match_indexes.Rdata')
save(regions, file = 'regions.Rdata')
save(coverageMatrix, file = 'coverageMatrix.Rdata')

## Session info:
options(width = 120)
devtools::session_info()
