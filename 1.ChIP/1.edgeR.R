library(liftOver)
library(openxlsx)
library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(cbmR)
library(SummarizedExperiment)
library(edgeR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library(org.Hs.eg.db)
require(ChIPseeker)

makeSumExp <- function(type)
{
  counts <- fread(paste0("rawData/", type, "_peakCounts.csv.gz"))
  binCounts <- fread(paste0("rawData/", type, "_genomicBins.csv.gz"))
  setcolorder(binCounts, colnames(counts))
  setnames(binCounts, str_remove(colnames(binCounts), paste0("_", type)))
  setnames(counts, str_remove(colnames(counts), paste0("_", type)))
  
  cDat <- data.table(id = colnames(counts))
  cDat[, c("Passage", "Treatment", "Replicate"):=tstrsplit(id, "_")]
  cDat[, Passage:=factor(Passage, levels = c("P5", "P6"))]
  cDat[, Treatment:=factor(Treatment, levels = c("Ctrl", "Palm", "TNFa"))]
  cDat[, Replicate:=factor(Replicate, levels = c("Rep1", "Rep2"))]
  
  setcolorder(counts, cDat$id)
  setcolorder(binCounts, cDat$id)
  
  annotation <- fread(paste0("rawData/", type, "_annotation.csv.gz"))
  rowRan <- with(annotation, 
                 GRanges(Chr, 
                         IRanges(Start, 
                                 End, 
                                 names = GeneID
                         ), 
                         seqinfo = seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)))
  
  sumExp <- SummarizedExperiment(
    assays = list(counts = as.matrix(counts)),
    rowRanges = rowRan,
    colData = cDat,
    metadata = list(binCounts = binCounts)
  )
  sumExp
}

k27ac <- makeSumExp("K27ac")
k4me1 <- makeSumExp("K4me1")

k27ac_gr <- rowRanges(k27ac)
k4me1_gr <- rowRanges(k4me1)

# Count overlaps
sum(overlapsAny(k27ac, k4me1))
length(k27ac)
# 77837 / 81074 K27ac peaks overlap a K4me1 peak

sum(overlapsAny(k4me1, k27ac))
length(k4me1)
# 51428 / 107405 K4me1 peaks overpal a K27ac peak

# We only care about situations where we have both K27ac and K4me1
k27ac <- subsetByOverlaps(k27ac, k4me1)
k4me1 <- subsetByOverlaps(k4me1, k27ac)

export(rowRanges(k27ac), con = "K27ac_overlapping_K4me1.bed")
export(rowRanges(k4me1), con = "K4me1_overlapping_K27ac.bed")

# Next remove active promoters
k4me3_hg19 <- import("h3k4me3_peaks/E121-H3K4me3.narrowPeak.gz")
# Somehow fdr ~ 0.46 is included in the dataset...
k4me3_hg19 <- k4me3_hg19[k4me3_hg19$qValue > -log10(0.05)]
ch <- import.chain("h3k4me3_peaks/hg19ToHg38.over.chain")
k4me3 <- liftOver(k4me3_hg19, ch) %>% unlist
length(k4me3)
# Far more than there are promoters. K4me3 can also be found in enhancers...
# Filter out only those in known promoters
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

annoDb <- "org.Hs.eg.db"

k4me3 <- k4me3 %>%
  unstrand() %>%
  keepStandardChromosomes() %>%
  annotatePeak(TxDb = txdb, annoDb = annoDb, tssRegion = c(3000, 1000)) %>%
  slot("anno")
k4me3 <- k4me3[k4me3$distanceToTSS > -3000 & k4me3$distanceToTSS < 1000]
activePromoters <- with(as.data.table(k4me3), GRanges(paste0("chr", geneChr), IRanges(geneStart, geneEnd), strand = ifelse(geneStrand == 1, "+", "-")))
seqlevels(activePromoters) <- str_replace(seqlevels(activePromoters), "chr23", "chrX")
seqlevels(activePromoters) <- str_replace(seqlevels(activePromoters), "chr24", "chrY")
activePromoters <- promoters(activePromoters, upstream = 3000, downstream = 1000)

k27ac <- subsetByOverlaps(k27ac, activePromoters, invert = TRUE)
export(rowRanges(k27ac), con = "enhancers.bed")

colScale <- scale_colour_manual(
  values = c(Ctrl = "#1b9e77",
             Palm = "#d95f02",
             TNFa = "#7570b3"),
  labels = c(Ctrl = "Control",
             Palm = "Palmitate",
             TNFa = "TNFa")
)

# Prepare the DGELists
source("../supportFunctions.R", local = TRUE)
dgeLists <- list(
  k27ac = chipPrepare(sumExp = k27ac, 
                      folder = "qcFigures/k27ac/", 
                      group = "Treatment", 
                      mdsLabel = "Passage",
                      mdsColour = "Treatment", 
                      colScale = colScale),
  
  # Create summary information on K4me1
  k4me1 = chipPrepare(sumExp = k4me1, 
                      folder = "qcFigures/k4me1/", 
                      group = "Treatment", 
                      mdsLabel = "Passage",
                      mdsColour = "Treatment", 
                      colScale = colScale)
)
# Run edgeR
edgerFn <- function(dgeList, folder)
{
  expSetup <- as.data.table(dgeList$samples)
  expSetup[, subgroup:=paste(Passage, Treatment, sep = "_")]
  design <- model.matrix(~ 0 + subgroup + Replicate, data = expSetup)
  colnames(design) <- str_replace_all(colnames(design), "Treatment|Passage|subgroup", "")
  colnames(design) <- str_replace(colnames(design), "Replicate", "")
  colnames(design) <- str_replace(colnames(design), ":", "_")
  
  y <- estimateDisp(dgeList, design = design, robust = TRUE)
  redirectToFile(y, file.path(folder, "BCV.pdf"), plotBCV)
  
  efit <- glmQLFit(y, design = design, robust = TRUE)
  redirectToFile(efit, file.path(folder, "QLDispersion.pdf"), plotQLDisp)
  
  colnames(design)
  ctrsts <- makeContrasts(
    Palmitate = (P5_Palm + P6_Palm)/2 - (P5_Ctrl + P6_Ctrl)/2,
    TNFa = (P5_TNFa + P6_TNFa)/2 - (P5_Ctrl + P6_Ctrl)/2,
    Replicate = Rep2,
    levels = design
  )
  
  res_ins <- apply(ctrsts, 2, . %>% 
                     glmQLFTest(glmfit = efit, contrast = .) %>% 
                     topTags(n = Inf, p.value = 1) %>% 
                     extract2("table") %>% 
                     as.data.table(keep.rownames = TRUE))
  res_ins
}

edgerRes <- list(
  k27ac = edgerFn(dgeLists$k27ac, "qcFigures/k27ac"),
  k4me1 = edgerFn(dgeLists$k4me1, "qcFigures/k4me1")
)

# Annotate ranges
annotateRanges <- function(resTable, txdb)
{
  resTable <- copy(resTable)
  annoDb <- "org.Hs.eg.db"
  
  ranges <- as(resTable, Class = "GRanges")
  setnames(resTable, "rn", "peakID")
  ranges <- ranges[order(ranges$PValue)]
  
  ranges %>%
    unstrand() %>%
    keepStandardChromosomes() %>%
    annotatePeak(TxDb = txdb, annoDb = annoDb) %>%
    slot("anno")
}

annotated <- list(
  k27ac = lapply(edgerRes$k27ac, annotateRanges, txdb = txdb),
  k4me1 = lapply(edgerRes$k4me1, annotateRanges, txdb = txdb)
)

write.xlsx(lapply(annotated$k27ac, as.data.frame), file = "enhancerPeaksAnnotated.xlsx", asTable = TRUE)
write.xlsx(lapply(annotated$k4me1, as.data.frame), file = "K4me1PeaksAnnotated.xlsx", asTable = TRUE)

annoDb <- "org.Hs.eg.db"
k27ac_gr %<>%
  unstrand() %>%
  keepStandardChromosomes(species = "Homo_sapiens") %>%
  annotatePeak(TxDb = txdb, annoDb = annoDb)
fwrite(k27ac_gr@annoStat, file = "k27ac_annotationStatistics.csv")

k4me1_gr %<>%
  unstrand() %>%
  keepStandardChromosomes(species = "Homo_sapiens") %>%
  annotatePeak(TxDb = txdb, annoDb = annoDb)
fwrite(k4me1_gr@annoStat, file = "k4me1_annotationStatistics.csv")
