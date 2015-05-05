####

library('derfinder')
library('GenomicRanges')
library('devtools')

### phenotype data
pd = read.delim("GD667.QCstats.masterfile.txt",	 as.is=TRUE)
pd = pd[,1:37]

# ## Load matching IDs
load("pMatch.Rdata")
pd$RailID = pMatch$railName[match(rownames(pd), pMatch$bgName)]

## Load regions data
load("/dcs01/ajaffe/Brain/derRuns/railDER/railGEU/fixSampleNames/regions.Rdata")
load("/dcs01/ajaffe/Brain/derRuns/railDER/railGEU/fixSampleNames/coverageMatrix.Rdata")

### filter to the same set of people ####
coverageMatrix = coverageMatrix[,pd$RailID]

#################
#### analysis ###

## expressed region 
load("/home/epi/ajaffe/GenomicStates/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")
gs = GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome
ensemblAnno = annotateRegions(regions,gs)
countTable = ensemblAnno$countTable

## annotation ####
dim(countTable)
annoClassList = list(strictExonic = 
	which(countTable[,"exon"] > 0 & countTable[,"intron"] == 0 &
		countTable[,"intergenic"] == 0),
	strictIntronic = 
	which(countTable[,"intron"] > 0 & countTable[,"exon"] == 0 &
		countTable[,"intergenic"] == 0),
	strictIntergenic = which(countTable[,"intergenic"] > 0 & countTable[,"exon"] == 0 &
    countTable[,"intron"] == 0),
	exonIntron = which(countTable[,"exon"] > 0 & countTable[,"intron"] > 0 &
		countTable[,"intergenic"] == 0))
sapply(annoClassList, length)
100*sapply(annoClassList, length)/nrow(countTable)

#################################
###### reproduce their figure ###

### joint modeling ####
y = log2(coverageMatrix+1)
colnames(y) = rownames(pd)
sumSqList = parallel::mclapply(1:nrow(y), function(i) {
	if(i %% 1000 == 0) cat(".")
	t(anova(lm(y[i,] ~ Population + RIN + RNAExtractionBatch + 
		RNAConcentration_ng.ul + RNAQuantityLibraryPrep_ng + 
		LibraryPrepDate + PrimerIndex + LibraryConcentrationMethod + 
		LibraryConcentration_ng.ul + BioanalyzerSize_bp + 
		LibraryQuantitySequencing_pM + ClusterKitBatch + 
		SequencingKitBatch + ClusterDensityPass + Lane, data=pd))[2])
},mc.cores=12)

ssOut = do.call("rbind", sumSqList)
rownames(ssOut) = NULL
bg = matrix(rep(rowSums(ssOut), ncol(ssOut)), 
	nc = ncol(ssOut),nrow = nrow(ssOut))
ssMat= ssOut / bg
lab = c("Population", "RIN value", "RNA extraction batch",
	"RNA concentration", "RNA quantity used",
	"Library preparation date", "Primer index",
	"Method concentration measure", "Library concentration",
	"Library size", "Library concentration used","Cluster kit",
	"Sequencing kit", "Cluster density", "Lane", "Residual variation")
save(ssMat, lab, file="ssMat_geuvadis.rda",compress=TRUE)

# load("ssMat_geuvadis.rda")

## overall boxplot
pdf("r2_boxplots_overall.pdf", h = 5, w = 12)
par(mar=c(9,5,2,2))
palette(brewer.pal(7, "Dark2"))
boxplot(100*ssMat,xaxt="n", ylim = c(0,90), 
	cex.axis=1.3,cex.lab=1.1, range=2,
	ylab="Percentage variance explained", cex=0.5)
text(1:ncol(ssMat)+0.2, y = -8, lab, xpd=TRUE, srt=45, pos=2)
text(x = 8.5, y= 80, "All Regions", cex=1.7)
for(i in seq(along=annoClassList)) {
	ii= annoClassList[[i]]
	boxplot(100*ssMat[ii,],xaxt="n", ylim = c(0,90), 
		cex.axis=1.3,cex.lab=1.1,range=2, col = i,
		ylab="Percentage variance explained", cex=0.5)
	text(1:ncol(ssMat)+0.1, y = -8, lab, xpd=TRUE, srt=45, pos=2)
	text(x = 8.5, y= 80, names(annoClassList)[i], cex=1.7)
}
dev.off()

## Reproducibility info
Sys.time() # date generated
options(width = 120)
session_info()
