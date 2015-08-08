####

library('derfinder')
library('GenomicRanges')
library('devtools')
library('RColorBrewer')
library('limma')

### phenotype data
pd = read.delim("GD667.QCstats.masterfile.txt",	 as.is=TRUE)
pd = pd[,1:37]

# ## Load matching IDs
load("/dcs01/ajaffe/Brain/derRuns/railDER/railGEU/fixSampleNames/pMatch.Rdata")
pd$RailID = pMatch$railName[match(rownames(pd), pMatch$bgName)]

## Load regions data
load("/dcs01/ajaffe/Brain/derRuns/railDER/resub/regionMatrix/regionMat-cut5.Rdata")
regions = unlist(GRangesList(lapply(regionMat, '[[', 'regions')))
coverageMatrix = do.call("rbind", lapply(regionMat, '[[', 'coverageMatrix'))
coverageMatrix = coverageMatrix[,pd$RailID] # put in order

#################
#### analysis ###

## expressed region 
load("/home/epi/ajaffe/GenomicStates/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")
gs = GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome
ensemblAnno = annotateRegions(regions,gs)
save(ensemblAnno, file = 'ensemblAnno.Rdata')
# load('ensemblAnno.Rdata')
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

quantile(width(regions))

sapply(annoClassList, function(ii) quantile(width(regions[ii])))

## Venn diagram: code modified from limma::vennDiagram
vennDiagram_custom <- function (object, include = "both", names = NULL, 
    mar = rep(1, 4), cex = c(1.5, 1, 0.7), lwd = 1, circle.col = NULL,
    counts.col = NULL, text.col = NULL, ...) 
{
    include <- as.character(include)
    LenInc <- min(length(include), 2)
    if (is(object, "VennCounts")) {
        include <- include[1]
        LenInc <- 1
    }
    else {
        if (LenInc > 1) 
            z2 <- vennCounts(object, include = include[2])[, 
                "Counts"]
        object <- vennCounts(object, include = include[1])
    }
    z <- object[, "Counts"]
    nsets <- ncol(object) - 1
    if (nsets > 5) 
        stop("Can't plot Venn diagram for more than 5 sets")
    VennZone <- object[, 1:nsets, drop = FALSE]
    VennZone <- apply(VennZone, 1, function(x) paste(x, sep = "", 
        collapse = ""))
    names(z) <- VennZone
    if (length(include) == 2) 
        names(z2) <- VennZone
    if (is.null(names)) 
        names <- colnames(object)[1:nsets]
    FILL.COL <- TRUE
    if (is.null(circle.col)) {
        circle.col <- par("col")
        FILL.COL <- FALSE
    }
    if (length(circle.col) < nsets) 
        circle.col <- rep(circle.col, length.out = nsets)
    if (is.null(counts.col)) 
        counts.col <- par("col")
    if (length(counts.col) < LenInc) 
        counts.col <- rep(counts.col, length.out = LenInc)
    if(is.null(text.col)) text.col <- rep('black', switch(nsets, counts.col[1], counts.col[1], 8))
    old.par <- par()$mar
    on.exit(par(mar = old.par))
    par(mar = mar)
    if (nsets <= 3) {
        plot(x = 0, y = 0, type = "n", xlim = c(-4, 4), ylim = c(-4, 
            4), xlab = "", ylab = "", axes = FALSE, ...)
            
        theta <- 2 * pi * (0:360)/360
        xcentres <- switch(nsets, 0, c(-1, 1), c(-1, 1, 0))
        ycentres <- switch(nsets, 0, c(0, 0), c(1, 1, -2)/sqrt(3))
        r <- 2
        xtext <- switch(nsets, -1.2, c(-1.2, 1.2), c(-1.2, 1.2, 
            0))
        ytext <- switch(nsets, 1.8, c(1.8, 1.8), c(3, 3, 
            -3.5))
        for (circle in 1:nsets) {
            if (!FILL.COL) 
                lines(xcentres[circle] + r * cos(theta), ycentres[circle] + 
                  r * sin(theta), lwd = lwd, col = circle.col[circle])
            if (FILL.COL) {
                RGB <- col2rgb(circle.col[circle])/255
                ALPHA <- 0.06
                RGB.ALP <- rgb(RGB[1, 1], RGB[2, 1], RGB[3, 1], 
                  alpha = ALPHA)
                polygon(xcentres[circle] + r * cos(theta), ycentres[circle] + 
                  r * sin(theta), border = circle.col[circle], 
                  lwd = lwd, col = RGB.ALP)
            }
            text(xtext[circle], ytext[circle], names[circle], 
                cex = cex * 1.3, col = circle.col[circle])
        }
        switch(nsets, rect(-3, -2.5, 3, 2.5), rect(-3, -2.5, 
            3, 2.5), rect(-3.9, -3.9, 3.9, 3.9))
        showCounts <- switch(nsets, function(counts, cex, adj, 
            col, leg) {
            text(2.3, -2.1, counts[1], cex = cex, col = col, 
                adj = adj)
            text(0, 0, counts[2], cex = cex, col = col, adj = adj)
        }, function(counts, cex, adj, col, leg) {
            text(2.3, -2.1, counts[1], cex = cex, col = col, 
                adj = adj)
            text(1.5, 0.1, counts[2], cex = cex, col = col, adj = adj)
            text(-1.5, 0.1, counts[3], cex = cex, col = col, 
                adj = adj)
            text(0, 0.1, counts[4], cex = cex, col = col, adj = adj)
        }, function(counts, cex, adj, col, leg) {
            text(3, -3, counts[1], cex = cex, col = col[1], adj = adj)
            text(0, -2.2, counts[2], cex = cex * 1.5, col = col[2], adj = adj)
            text(2, 1, counts[3], cex = cex * 1.5, col = col[3], adj = adj)
            text(1.3, -0.5, counts[4], cex = cex, col = col[4], 
                adj = adj)
            text(-2, 1, counts[5], cex = cex * 1.5, col = col[5], adj = adj)
            text(-1.3, -0.5, counts[6], cex = cex * 1.3, col = col[6], 
                adj = adj)
            text(0, 1.3, counts[7], cex = cex, col = col[7], adj = adj)
            text(0, 0, counts[8], cex = cex, col = col[8], adj = adj)
        })
        if (LenInc == 1) 
            adj <- c(0.5, 0.5)
        else adj <- c(0.5, 0)
        showCounts(counts = z, cex = cex[1], adj = adj, col = text.col, 
            leg = include[1])
        return(invisible())
    }
}



pdf('ensemblVenn.pdf', width = 10, height = 10)
venn_col <- brewer.pal(7, "Dark2")[c(1, 3, 2, 4)]
vennDiagram_custom(vennCounts(countTable > 0), 
    main = 'Expressed Regions overlap with Ensembl v75 features', cex.main = 2,
    circle.col = venn_col[1:3], lwd = 1.5, cex = 2, mar = c(0, 0, 2, 0),
    text.col = c('black', venn_col[3:2], 'black', venn_col[c(1, 4)], 'black',
        'black')#, oma = rep(0, 4), pty = 'm'
)
dev.off()

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

## Without population term
sumSqList2 = parallel::mclapply(1:nrow(y), function(i) {
	if(i %% 1000 == 0) cat(".")
	t(anova(lm(y[i,] ~ RIN + RNAExtractionBatch + 
		RNAConcentration_ng.ul + RNAQuantityLibraryPrep_ng + 
		LibraryPrepDate + PrimerIndex + LibraryConcentrationMethod + 
		LibraryConcentration_ng.ul + BioanalyzerSize_bp + 
		LibraryQuantitySequencing_pM + ClusterKitBatch + 
		SequencingKitBatch + ClusterDensityPass + Lane, data=pd))[2])
}, mc.cores=12)

ssOut2 = do.call("rbind", sumSqList2)
rownames(ssOut2) = NULL
bg2 = matrix(rep(rowSums(ssOut2), ncol(ssOut2)), 
	nc = ncol(ssOut2),nrow = nrow(ssOut2))
ssMat2 = ssOut2 / bg2
lab2 = c("RIN value", "RNA extraction batch",
	"RNA concentration", "RNA quantity used",
	"Library preparation date", "Primer index",
	"Method concentration measure", "Library concentration",
	"Library size", "Library concentration used","Cluster kit",
	"Sequencing kit", "Cluster density", "Lane", "Residual variation")
save(ssMat2, lab2, file="ssMat_geuvadis_noPop.rda", compress=TRUE)

## overall boxplot no pop
pdf("r2_boxplots_overall_noPop.pdf", h = 5, w = 12)
par(mar=c(9,5,2,2))
palette(brewer.pal(7, "Dark2"))
boxplot(100*ssMat2,xaxt="n", ylim = c(0,90), 
	cex.axis=1.3,cex.lab=1.1, range=2,
	ylab="Percentage variance explained", cex=0.5)
text(1:ncol(ssMat2)+0.2, y = -8, lab2, xpd=TRUE, srt=45, pos=2)
text(x = 8.5, y= 80, "All Regions", cex=1.7)
for(i in seq(along=annoClassList)) {
	ii= annoClassList[[i]]
	boxplot(100*ssMat2[ii,],xaxt="n", ylim = c(0,90), 
		cex.axis=1.3,cex.lab=1.1,range=2, col = i,
		ylab="Percentage variance explained", cex=0.5)
	text(1:ncol(ssMat2)+0.1, y = -8, lab2, xpd=TRUE, srt=45, pos=2)
	text(x = 8.5, y= 80, names(annoClassList)[i], cex=1.7)
}
dev.off()

## Reproducibility info
Sys.time() # date generated
options(width = 120)
session_info()
