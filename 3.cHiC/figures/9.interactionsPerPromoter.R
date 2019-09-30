library(stringr)
library(data.table)
library(ggplot2)
library(magrittr)
library(openxlsx)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(InteractionSet)
load("../selectedInteractions.RData")

promoters <- import("../promoters_hg38_fixed.bed.gz")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

annotateFn <- function(ranges, txdb)
{
  annoDb <- "org.Hs.eg.db"
  
  ranges %>%
    unstrand() %>%
    keepStandardChromosomes(species = "Homo_sapiens") %>%
    annotatePeak(TxDb = txdb, annoDb = annoDb, addFlankGeneInfo = TRUE, flankDistance = 10000) %>%
    slot("anno")
}

sortEnds <- function(x, isSecond = 0)
{
  lapply(anchors(x), subset, is.second == isSecond) %>% 
    GRangesList %>% 
    unlist %>%
    unname
}

genes <- sortEnds(viewpoints) %>% annotateFn(txdb = txdb)
enhancerEnds <- sortEnds(viewpoints, isSecond = 1)

#### Convert the geneIDs to ENSEMBL IDs
txIdToENSEMBL <- org.Hs.egENSEMBL %>%
  {.[mappedkeys(.)]} %>%
  as.list %>%
  unlist2 %>%
  {data.table(tx_id = names(.), ENSEMBL = .)} %T>%
  setkey(tx_id)

txIdToSYMBOL <- org.Hs.egSYMBOL %>%
{.[mappedkeys(.)]} %>%
  as.list %>%
  unlist2 %>%
  {data.table(tx_id = names(.), SYMBOL = .)} %T>%
  setkey(tx_id)

ENSEMBLtoSYMBOL <- txIdToENSEMBL[txIdToSYMBOL, nomatch = FALSE][, list(SYMBOL), key = "ENSEMBL"]

# Get the new list of nearby genes (for every gene remove situations where the flanking gene is the same as the main annotation)
genes_close <- subset(genes, distanceToTSS < 3000)
flankGenes <- genes_close$flank_geneIds %>% 
  str_split(pattern = ";") %>%
  lapply(unique)

nClose <- lapply(flankGenes, length) %>% unlist
table(nClose) %>% plot
unlist(flankGenes) %>% table %>% table %>% sort %>% plot(xlim = c(0, 20))
genes_close$geneId %>% table %>% table %>% sort %>% plot(xlim = c(0, 20))

updatedNearbyGenes <- Map(union, as.list(genes_close$geneId), flankGenes) %>%
  unlist

updatedNearbyGenes <- genes_close$geneId

genesByNumberOfInteractions <- txIdToENSEMBL[updatedNearbyGenes, ENSEMBL, mult = "first"] %>%
  na.omit %>% 
  table %>% 
  {.[order(., decreasing = TRUE)]} %>%
  as.data.table %>%
  setnames(c("ENSEMBL", "N"))

genesByNumberOfInteractions[, SYMBOL:=ENSEMBLtoSYMBOL[ENSEMBL, SYMBOL, mult = "first"]]

write.csv(genesByNumberOfInteractions, file = "genesByNumberOfInteractions.csv", row.names = FALSE)

dens <- density(genesByNumberOfInteractions$N)

modeInteractions <- dens$x[which.max(dens$y)]
medianInteractions <- median(genesByNumberOfInteractions$N, na.rm = TRUE)
meanInteractions <- mean(genesByNumberOfInteractions$N, na.rm = TRUE)


qplot(as.vector(genesByNumberOfInteractions$N), geom="histogram", bins = 50) +
  scale_x_continuous(name = "Number of interactions", limits = c(0, 50)) +
  scale_y_continuous(name = element_blank()) + 
  theme_bw()

ggsave(filename = "interactionsPerPromoter.pdf")


### To calculate the number of interaction pr. enhancer simply count how many overlap the fragments
enhancers <- import("../../1.ChIP/enhancers.bed")
qplot(countOverlaps(enhancers, enhancerEnds), geom = "histogram", bins = 50) +
  scale_x_continuous(name = "Number of interactions", limits = c(-1, 50)) +
  scale_y_continuous(name = element_blank()) + 
  theme_cbmr()

ggsave("Fig3B_enhancers.pdf")

