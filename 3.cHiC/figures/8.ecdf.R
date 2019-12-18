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
rnaSeq <- readXlsxFile("../../2.RNA/edgeR_results/RNA.xlsx")

enhancer <- readXlsxFile("../../1.ChIP/enhancerPeaksAnnotated.xlsx")
lapply(enhancer, setkey, "rn") %>% invisible()
enhancers_signif <- list(
  Palmitate = enhancer$Palmitate[FDR < 0.01, ],
  TNFa = enhancer$TNFa[FDR < 0.01, ]
)

load("../lookupTables/combinedResults.RData")
p2e <- promoters2enhancers
enhancers <- enhancers_signif$Palmitate
rna <- rnaSeq$Palmitate

preparePlotData <- function(enhancers, rna, p2e)
{
  p2e <- copy(p2e)
  p2e <- p2e[, .(rn = unlist(enhancerID)), by = "promoterIDs"]
  cmb <- enhancers[p2e, nomatch = FALSE, on = "rn"]
  helper <- function(x)
  {
    x <- sign(x)
    if (all(x == 1))
    {
      return("Up")
    } else if (all(x == -1))
    {
      return("Down")
    } else {
      return("Both")
    }
    stop("How did I get here?", x)
  }
  directions <- cmb[, helper(logFC), by = "promoterIDs"]
  
  pD <- copy(rna)
  pD[, group:="None"]
  setkey(pD, "ENSEMBL")
  pD[directions$promoterIDs, group:=directions$V1]
  pD
}
pD <- list(
  Palmitate = preparePlotData(enhancers_signif$Palmitate, rnaSeq$Palmitate, p2e = promoters2enhancers),
  TNFa = preparePlotData(enhancers_signif$TNFa, rnaSeq$TNFa, p2e = promoters2enhancers)
)


colScale <- scale_colour_manual(name = element_blank(), values = c("None" = "#4daf4a", "Up" = "#e41a1c", "Down" = "#377eb8", "Both" = "#984ea3"))
baseplot <- ggplot(mapping = aes(x = logFC, colour = group)) + 
  stat_ecdf(size = 2) +
  theme_bw(base_size = 20) +
  geom_vline(xintercept = 0, lty = 2) +
  coord_cartesian(xlim = c(-1, 1)) +
  xlab("Estimated logFC") +
  scale_y_continuous(name = "Cumulative percent", labels = scales::percent) +
  colScale +
  theme(legend.justification = c(1, 0), legend.position = c(0.9, 0.1))

width = 20
height = 10

tiff("ecdfPlots_withHiC/TNFa_all.tiff", width = width, height = height, units = "cm", compression = "lzw+p", type = "cairo", res = 1200)
p <- baseplot %+% pD$TNFa[group %in% c("None", "Both", "Down", "Up")] +
  coord_cartesian(xlim = c(-1, 1)) 
print(p)
dev.off()

tiff("ecdfPlots_withHiC/TNFa_up_down.tiff", width = width, height = height, units = "cm", compression = "lzw+p", type = "cairo", res = 1200)
p <- baseplot %+% pD$TNFa[group %in% c("None", "Down", "Up")] +
  coord_cartesian(xlim = c(-1, 1)) 
print(p)
dev.off()

tiff("ecdfPlots_withHiC/TNFa_up.tiff", width = width, height = height, units = "cm", compression = "lzw+p", type = "cairo", res = 1200)
p <- baseplot %+% pD$TNFa[group %in% c("None", "Up")] +
  coord_cartesian(xlim = c(-1, 1)) 
print(p)
dev.off()

tiff("ecdfPlots_withHiC/TNFa_down.tiff", width = width, height = height, units = "cm", compression = "lzw+p", type = "cairo", res = 1200)
p <- baseplot %+% pD$TNFa[group %in% c("None", "Down")] +
  coord_cartesian(xlim = c(-1, 1)) 
print(p)
dev.off()

tiff("ecdfPlots_withHiC/TNFa_both.tiff", width = width, height = height, units = "cm", compression = "lzw+p", type = "cairo", res = 1200)
p <- baseplot %+% pD$TNFa[group %in% c("None", "Both")] +
  coord_cartesian(xlim = c(-1, 1)) 
print(p)
dev.off()



tiff("ecdfPlots_withHiC/Palmitate_all.tiff", width = width, height = height, units = "cm", compression = "lzw+p", type = "cairo", res = 1200)
p <- baseplot %+% pD$Palmitate[group %in% c("None", "Both", "Down", "Up")]
print(p)
dev.off()

tiff("ecdfPlots_withHiC/Palmitate_up_down.tiff", width = width, height = height, units = "cm", compression = "lzw+p", type = "cairo", res = 1200)
p <- baseplot %+% pD$Palmitate[group %in% c("None", "Down", "Up")]
print(p)
dev.off()

tiff("ecdfPlots_withHiC/Palmitate_up.tiff", width = width, height = height, units = "cm", compression = "lzw+p", type = "cairo", res = 1200)
p <- baseplot %+% pD$Palmitate[group %in% c("None", "Up")]
print(p)
dev.off()

tiff("ecdfPlots_withHiC/Palmitate_down.tiff", width = width, height = height, units = "cm", compression = "lzw+p", type = "cairo", res = 1200)
p <- baseplot %+% pD$Palmitate[group %in% c("None", "Down")]
print(p)
dev.off()

tiff("ecdfPlots_withHiC/Palmitate_both.tiff", width = width, height = height, units = "cm", compression = "lzw+p", type = "cairo", res = 1200)
p <- baseplot %+% pD$Palmitate[group %in% c("None", "Both")]
print(p)
dev.off()

##### Get P-values
signifTester <- function(dt, group)
{
  ks.test(dt[group, logFC, on = "group"], dt["None", logFC, on = "group"])
}
tests <- c("Up", "Down", "Both")
names(tests) <- tests
lapply(tests, signifTester, dt = pD$Palmitate)

lapply(tests, signifTester, dt = pD$TNFa)
