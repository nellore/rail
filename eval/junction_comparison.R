###
library(GenomicRanges)
splitit = function(x) split(seq(along=x),x) # splits into list

##### make tophat count table
source("junctionCount.R")
load("/home/epi/ajaffe/Lieber/Projects/RNAseq/theControlEqtl/rdas/geuvadis_junctions.rda")
juncCounts = junctionCount(theJunctions, maxCores=12)

## read in ensembl Junctions
load("ensembl_v75_junction_annotation.rda")

### read in tophat jxns
jMapTop = juncCounts$anno
jCountTop = as.data.frame(juncCounts$countDF)
jCountTop = as.matrix(jCountTop)

jMapTop$totalReads = rowSums(jCountTop)

jMapTop$inEnsembl = countOverlaps(jMapTop, theJunctions, type="equal") > 0
jMapTop$inEnsemblStart = countOverlaps(jMapTop, theJunctions, type="start") > 0
jMapTop$inEnsemblEnd = countOverlaps(jMapTop, theJunctions, type="end") > 0
table(jMapTop$inEnsembl)
jMapTop$meanReads = jMapTop$totalReads / ncol(jCountTop)
jMapTop$code = ifelse(jMapTop$inEnsembl, "InEns", 
	ifelse(jMapTop$inEnsemblStart & jMapTop$inEnsemblEnd, "NovelTrans",
	ifelse(jMapTop$inEnsemblStart | jMapTop$inEnsemblEnd, "NovelJxn", "NewJxn")))

## read in Rail jxns
railJxn = read.table("/dcs01/ajaffe/GEUVADIS/intron_table.tsv",
	sep="\t",header=FALSE)
railJxnMat =  read.table("/dcs01/ajaffe/GEUVADIS/GEUVADIS_introns_v4.04.2015.tsv.gz",
	header=TRUE,sep="\t",row.names=1)

jMapRail = with(railJxn, GRanges(V1, IRanges(V2, V3-1)))
names(jMapRail) = paste0(seqnames(jMapRail), ":", 
	start(jMapRail), "-", end(jMapRail), "(*)")
jMapRail$totalReads = rowSums(railJxnMat)
jMapRail$inEnsembl = countOverlaps(jMapRail, theJunctions, type="equal") > 0
jMapRail$inEnsemblStart = countOverlaps(jMapRail, theJunctions, type="start") > 0
jMapRail$inEnsemblEnd = countOverlaps(jMapRail, theJunctions, type="end") > 0
jMapRail$meanReads = jMapRail$totalReads / ncol(jCountTop)
jMapRail$code = ifelse(jMapRail$inEnsembl, "InEns", 
	ifelse(jMapRail$inEnsemblStart & jMapRail$inEnsemblEnd, "NovelTrans",
	ifelse(jMapRail$inEnsemblStart | jMapRail$inEnsemblEnd, "NovelJxn", "NewJxn")))

	
###### compare #####	
length(jMapTop)	# number of jxns
length(jMapRail) # number of jxns	

topFilter = rowMeans(jCountTop > 0) > 0.05 | 
	rowSums(jCountTop > 4) > 0
mean(topFilter)

######## novel transcripts #######
## gene ranges ###################
jList = split(theJunctions, theJunctions$ensemblID)
geneRanges = unlist(range(jList))
geneRanges = geneRanges[order(geneRanges)]
geneRanges$coordRange = paste0(seqnames(geneRanges), ":",
	start(geneRanges), "-", end(geneRanges))
hla = GRanges("chr6", IRanges(29691116,33054976)) ## hla region
geneRanges$inHLA = countOverlaps(geneRanges, hla) > 0

## check rail
jMapRailHi = jMapRail[jMapRail$meanReads > 10]
oo = findOverlaps(jMapRailHi, geneRanges)
jRailList = CharacterList(split(jMapRailHi$code[queryHits(oo)], 
	names(geneRanges)[subjectHits(oo)]))
jRailTab = as.data.frame(sapply(unique(jMapRailHi$code), 
	function(x) sum(jRailList == x)))
jRailTab$Symbol = theJunctions$symbol[
	match(rownames(jRailTab), theJunctions$ensemblID)]
jRailTab$Coords = geneRanges$coordRange[match(rownames(jRailTab), names(geneRanges))]
jRailTab$inHLA = geneRanges$inHLA[match(rownames(jRailTab), names(geneRanges))]

## pseudogene
pseudo = read.delim("/nexsan2/disk3/ajaffe/RNASeq/PseudoPipeHuman61.txt",
	as.is=TRUE)[,-22]
jRailTab$hasPseudo = rownames(jRailTab) %in% pseudo$Parent.Gene
write.csv(jRailTab, file="junctionBreakdown_byLocusRange.csv")
	
#### any genes with many novel transcripts?
jRailTab2 = jRailTab[order(jRailTab$NovelTrans,decreasing=TRUE),]	
mean(jRailTab2$hasPseudo[1:50])
# jRailTab2 = jRailTab2[!is.na(jRailTab2$Symbol),]
mean(jRailTab2$inHLA[1:50])
hlaRailIndex = grep("HLA", jRailTab2$Symbol)
jRailTab2[hlaRailIndex,]

## tophat results	
jMapTopHi = jMapTop[topFilter & jMapTop$meanReads > 10]
oo = findOverlaps(jMapTopHi, geneRanges)
jTopList = CharacterList(split(jMapTopHi$code[queryHits(oo)], 
	names(geneRanges)[subjectHits(oo)]))
jTopTab = as.data.frame(sapply(unique(jMapRailHi$code), 
	function(x) sum(jTopList == x)))
jTopTab$Symbol = theJunctions$symbol[
	match(rownames(jTopTab), theJunctions$ensemblID)]
jTopTab$Coords = geneRanges$coordRange[match(rownames(jTopTab), names(geneRanges))]
jTopTab$inHLA = geneRanges$inHLA[match(rownames(jTopTab), names(geneRanges))]
jTopTab$hasPseudo = rownames(jTopTab) %in% pseudo$Parent.Gene

#### any genes with many novel transcripts?
jTopTab2 = jTopTab[order(jTopTab$NovelTrans,decreasing=TRUE),]	
mean(jTopTab2$inHLA[1:50])
jTopTab2 = jTopTab2[!is.na(jTopTab2$Symbol),]
head(jTopTab2, 25)
hlaTopIndex =  grep("HLA", jTopTab2$Symbol)


############################
### proportion table #######
tabOut = rbind(table(jMapTop$code),
	table(jMapTop$code[topFilter]),
	table(jMapRail$code))
rownames(tabOut) = c("TopHat 2", "TopHat 2 Filtered", "Rail-RNA")
matOut = matrix(paste0(tabOut, " (", 100*signif(prop.table(tabOut, 1),3),"%)"),
	nr = nrow(tabOut), nc = ncol(tabOut), dimnames=dimnames(tabOut))
matOut = cbind(matOut, rowSums(tabOut))
colnames(matOut)[5] = "Total"
write.csv(matOut, file="junctionTable_countsAndPerc.csv")

### do venn diagram
library(limma)
uJxn = c(jMapRail, jMapTop)
uJxn = uJxn[!duplicated(uJxn)]
mcols(uJxn) = mcols(uJxn)[-c(1,5)]
vennTab = matrix(FALSE, nr = length(uJxn), nc = 3)
rownames(vennTab) = paste0(seqnames(uJxn), ":", 
	start(uJxn), "-", end(uJxn), "(*)")
colnames(vennTab) = c("TopHat 2", "TopHat 2 Filtered", "Rail-RNA")
vennTab[rownames(vennTab) %in% names(jMapTop), 1] = TRUE
vennTab[rownames(vennTab) %in% names(jMapTop[topFilter]), 2] = TRUE
vennTab[rownames(vennTab) %in% names(jMapRail), 3] = TRUE

pdf("vennDiagrams_junctionClasses.pdf",h=6,w=7)
vennDiagram(vennCounts(vennTab),mar=rep(0.1, 4), cex=1.1)
mtext("All Classes", line=1,cex=2)
vennDiagram(vennCounts(vennTab[uJxn$code == "InEns",]),
	mar=rep(0.1, 4),cex=1.1)
mtext("In Ensembl", line=1,cex=2)
vennDiagram(vennCounts(vennTab[uJxn$code == "NovelTrans",]),
	mar=rep(0.1, 4), cex=1.1)
mtext("Novel Transcript", line=1,cex=2)
vennDiagram(vennCounts(vennTab[uJxn$code == "NovelJxn",]),
	mar=rep(0.1, 4),cex=1.1)
mtext("Novel Junction", line=1,cex=2)
vennDiagram(vennCounts(vennTab[uJxn$code == "NewJxn",]),
	mar=rep(0.1, 4),cex=1.1)
mtext("New Junction", line=1,cex=2)
dev.off()

###############################
####  make boxplots ######

pdf("boxplots_totalMapped_byClass.pdf",h=5,w=8)
par(mar=c(5,6,2,2))
palette(RColorBrewer::brewer.pal(8,"Set1"))
boxplot(log2(jMapTop$totalReads+1) ~ jMapTop$code,range=2,
	col=1:4, cex.axis=1.6,cex.lab=1.6, cex.main=1.6,
	main = "TopHat 2",ylab="log2(Total Reads+1)",ylim=c(0,25))
boxplot(log2(jMapTop$totalReads+1) ~ jMapTop$code,range=2,
	col=1:4, cex.axis=1.6,cex.lab=1.6, cex.main=1.6,
	subset= topFilter, main = "TopHat 2 Filtered",
	ylab="log2(Total Reads+1)",ylim=c(0,25))
boxplot(log2(jMapRail$totalReads+1) ~ jMapRail$code,range=2,
	col=1:4, cex.axis=1.6,cex.lab=1.6, cex.main=1.6,
	subset= topFilter, main = "Rail-RNA",
	ylab="log2(Total Reads+1)",ylim=c(0,25))
dev.off()