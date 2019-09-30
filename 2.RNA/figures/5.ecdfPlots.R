library(biomaRt)
library(data.table)
library(magrittr)
library(ggplot2)
library(ChIPseeker)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(stringr)
library(clusterProfiler)
library(openxlsx)

readXlsxFile <- function(file) 
{
  x <- getSheetNames(file)
  names(x) <- x
  lapply(x, . %>% read.xlsx(xlsxFile = file) %>% data.table)
}
rnaSeq <- readXlsxFile("../edgeR_results/RNA.xlsx")

ensembl <- useMart('ENSEMBL_MART_ENSEMBL', 
                   dataset = 'hsapiens_gene_ensembl') %>%
  getBM(attributes= c("ensembl_gene_id",
                      "chromosome_name",
                      "start_position",
                      "strand",
                      "end_position"
  ), 
  mart = .) %>% 
  data.table
ensembl[, chromosome_name:=paste0("chr", chromosome_name)]
setnames(ensembl, c("ENSEMBL", "seqnames", "start", "strand", "end"))
ensembl <- ensembl[ENSEMBL %in% rnaSeq$Palmitate$ENSEMBL]
ensembl[, strand:=ifelse(strand == "-1", "-", "+")]
geneTss <- as(ensembl, "GRanges") %>% 
  keepStandardChromosomes(species = "Homo_sapiens", pruning.mode = "coarse") %>%
  promoters(upstream = 0, downstream = 0)

enhancer <- readXlsxFile("../../1.ChIP/enhancerPeaksAnnotated.xlsx")
enhancers_signif <- list(
  Palmitate = enhancer$Palmitate[FDR < 0.01, ] %>% as("GRanges"),
  TNFa = enhancer$TNFa[FDR < 0.01, ] %>% as("GRanges")
)


preparePlotData <- function(enhancers, geneTss, rnaSeq, maxGap)
{
  enhancerUp <- findOverlaps(geneTss, enhancers[enhancers$logFC > 0], maxgap = maxGap)
  enhancerDown <- findOverlaps(geneTss, enhancers[enhancers$logFC < 0], maxgap = maxGap)
  both <- intersect(queryHits(enhancerUp), queryHits(enhancerDown))
  
  pD <- rnaSeq[geneTss$ENSEMBL, on = "ENSEMBL"]
  pD[, group:="None"]
  pD[geneTss[queryHits(enhancerDown)]$ENSEMBL, group:="Down", on = "ENSEMBL"]
  pD[geneTss[queryHits(enhancerUp)]$ENSEMBL, group:="Up", on = "ENSEMBL"]
  pD[geneTss[both]$ENSEMBL, group:="Both", on = "ENSEMBL"]
  pD
}
pD <- list(
  Palmitate = preparePlotData(enhancers_signif$Palmitate, geneTss, rnaSeq$Palmitate, maxGap = 50000),
  TNFa = preparePlotData(enhancers_signif$TNFa, geneTss, rnaSeq$TNFa, maxGap = 50000)
)
    

colScale <- scale_colour_manual(name = element_blank(), values = c("None" = "#4daf4a", "Up" = "#e41a1c", "Down" = "#377eb8", "Both" = "#984ea3"))
baseplot <- ggplot(mapping = aes(x = logFC, colour = group)) + 
  stat_ecdf() +
  theme_bw() +
  geom_vline(xintercept = 0, lty = 2) +
  coord_cartesian(xlim = c(-1, 1)) +
  xlab("Estimated logFC") +
  scale_y_continuous(name = "Cumulative percent", labels = scales::percent) +
  colScale

tiff("ecdfPlots_noHiC/TNFa_all.tiff", width = 20, height = 20, units = "cm", compression = "lzw+p", type = "cairo", res = 1200)
p <- baseplot %+% pD$TNFa[group %in% c("None", "Both", "Down", "Up")] +
  coord_cartesian(xlim = c(-1.5, 3)) 
print(p)
dev.off()

tiff("ecdfPlots_noHiC/TNFa_up_down.tiff", width = 20, height = 20, units = "cm", compression = "lzw+p", type = "cairo", res = 1200)
p <- baseplot %+% pD$TNFa[group %in% c("None", "Down", "Up")] +
  coord_cartesian(xlim = c(-1.5, 3)) 
print(p)
dev.off()

tiff("ecdfPlots_noHiC/TNFa_up.tiff", width = 20, height = 20, units = "cm", compression = "lzw+p", type = "cairo", res = 1200)
p <- baseplot %+% pD$TNFa[group %in% c("None", "Up")] +
  coord_cartesian(xlim = c(-1.5, 3)) 
print(p)
dev.off()

tiff("ecdfPlots_noHiC/TNFa_down.tiff", width = 20, height = 20, units = "cm", compression = "lzw+p", type = "cairo", res = 1200)
p <- baseplot %+% pD$TNFa[group %in% c("None", "Down")] +
  coord_cartesian(xlim = c(-1.5, 3)) 
print(p)
dev.off()



tiff("ecdfPlots_noHiC/Palmitate_all.tiff", width = 20, height = 20, units = "cm", compression = "lzw+p", type = "cairo", res = 1200)
p <- baseplot %+% pD$Palmitate[group %in% c("None", "Both", "Down", "Up")] 
print(p)
dev.off()

tiff("ecdfPlots_noHiC/Palmitate_up_down.tiff", width = 20, height = 20, units = "cm", compression = "lzw+p", type = "cairo", res = 1200)
p <- baseplot %+% pD$Palmitate[group %in% c("None", "Down", "Up")] 
print(p)
dev.off()

tiff("ecdfPlots_noHiC/Palmitate_up.tiff", width = 20, height = 20, units = "cm", compression = "lzw+p", type = "cairo", res = 1200)
p <- baseplot %+% pD$Palmitate[group %in% c("None", "Up")] 
print(p)
dev.off()

tiff("ecdfPlots_noHiC/Palmitate_down.tiff", width = 20, height = 20, units = "cm", compression = "lzw+p", type = "cairo", res = 1200)
p <- baseplot %+% pD$Palmitate[group %in% c("None", "Down")] 
print(p)
dev.off()
