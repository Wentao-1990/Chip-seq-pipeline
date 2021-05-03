## loading packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)

files <- list(o3hoxa9_1="./o3hoxa9_rep1_peaks.narrowPeak", o3hoxa9_2="./o3hoxa9_rep2_peaks.narrowPeak", o3hoxa9_3="./o3hoxa9_rep3_peaks.narrowPeak")

peak <- readPeakFile(files[[1]])

##Chip peaks coverage plot##
pdf(file = "covplotforallchr.pdf",
    width = 16, height = 9,  
    bg = "transparent")
#covplot(peak, weightCol="V5")
covplot(peak, weightCol="V5",chrs=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
   "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"))
dev.off()

#pdf(file = "covplotforallchr.pdf",
#    width = 16, height = 9,  
#    bg = "transparent")
#covplot(peak, weightCol="V5", chrs=c("chr17", "chr18"), xlim=c(4.5e7, 5e7))

##Profile of ChIP peaks binding to TSS regions##
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)

##heatmap of ChIP binding to TSS region##
pdf(file = "heatmapforonetss.pdf",
    width = 16, height = 9,  
    bg = "transparent")
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
dev.off()

##a one step function to generate this figure###
#peakHeatmap(files[[4]], TxDb=txdb, upstream=3000, downstream=3000, color="red")

##Average Profile of Chip peaks binding to TSS region##
pdf(file = "profierplotforonetss.pdf",
    width = 16, height = 9,  
    bg = "transparent")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()

##a one step function from bed file to average profile plot##
#plotAvgProf2(files[[4]], TxDb=txdb, upstream=3000, downstream=3000,
#             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

###confidence interval estimated by bootstrap method##
pdf(file = "profierplotforonetsswithcon.pdf",
    width = 16, height = 9,  
    bg = "transparent")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)
dev.off()

###Peak Annotation##
peakAnno <- annotatePeak(files[[1]], tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

##Visualize Genomic Annotation##
pdf(file = "peakannpieplotforone.pdf",
    width = 16, height = 9,  
    bg = "transparent")
plotAnnoPie(peakAnno)
dev.off()
##Visuliza distribution of TF-binding loci relative to TSS##
pdf(file = "peakdisplotforone.pdf",
    width = 16, height = 9,  
    bg = "transparent")
plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")
dev.off()

##Functional enrichment analysis##
library(ReactomePA)

pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId)

gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene)

pdf(file = "reactenrichforone.pdf",
    width = 16, height = 9,  
    bg = "transparent")
dotplot(pathway2)
dev.off()

##Chip peak data set comparison###
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
pdf(file = "profierplotforgroup.pdf",
    width = 16, height = 9,  
    bg = "transparent")
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
dev.off()
##peak heatmaps##
pdf(file = "pheatmapplotforgroup.pdf",
    width = 16, height = 9,  
    bg = "transparent")
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
dev.off()

##Chip peak annotation comparison##

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)
pdf(file = "peakannobarforgroup.pdf",
    width = 16, height = 9,  
    bg = "transparent")
plotAnnoBar(peakAnnoList)
dev.off()
pdf(file = "peakannotssforgroup.pdf",
    width = 16, height = 9,  
    bg = "transparent")
plotDistToTSS(peakAnnoList)
dev.off()


##outputAnnoList file###
lapply(peakAnnoList, function(i) write.table(i@anno,file=paste0(i,'.txt',sep=""),sep="\t"))

##Functional profiles comparison##

genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
pdf(file = "enrichplotforgroup.pdf",
    width = 16, height = 9,  
    bg = "transparent")
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
dev.off()




##Overlap of peaks and annoated genes##

genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
pdf(file = "overlapplotforgroup.pdf",
    width = 16, height = 9,  
    bg = "transparent")
vennplot(genes)
dev.off()

##Statistical testing of ChIP seq overlap##
#Shuffle genome coordinatiion##
#p <- GRanges(seqnames=c("chr1", "chr3"),
#             ranges=IRanges(start=c(1, 100), end=c(50, 130)))
#shuffle(p, TxDb=txdb)

##Peak overlap enrichment analysis##
#enrichPeakOverlap(queryPeak     = files[[5]],
#                  targetPeak    = unlist(files[1:4]),
#                  TxDb          = txdb,
#                  pAdjustMethod = "BH",
#                  nShuffle      = 50,
#                  chainFile     = NULL,
#                  verbose       = FALSE)
