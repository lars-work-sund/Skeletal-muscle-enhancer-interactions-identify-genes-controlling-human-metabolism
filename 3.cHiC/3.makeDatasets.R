options(java.parameters = "-Xmx4096m")  ## memory set to 4 GB
library(diffHic)
require(BSgenome.Hsapiens.UCSC.hg38)
require(edgeR)
library(csaw)
library(magrittr)
library(data.table)
library(openxlsx)
library(stringr)
### This section requires the raw data which takes up many gigabytes of space.
### To save space only the final result is included in the data uploaded.


# hs.frag <- cutGenome(BSgenome.Hsapiens.UCSC.hg38, "AAGCTT", 4)
# hs.param <- pairParam(hs.frag)
# ## ---- Read files and reformat
# fls <- list.files("bamFiles/", ".bam", full.names = TRUE)
# names(fls) <- basename(fls)
# diagnostics <- Map(preparePairs, 
#                    bam = fls, 
#                    file = sub(".bam", ".h5", fls),
#                    MoreArgs = list(param = hs.param)
#                    )
# 
# min.inward <- 1000
# min.outward <- 25000
# 
# pruning <- Map(prunePairs,
#     file.in = sub(".bam", ".h5", fls),
#     file.out = sub(".bam", "_trimmed.h5", fls),
#     MoreArgs = list(param = hs.param,
#                     min.inward = min.inward,
#                     min.outward = min.outward,
#                     max.frag=600))

## ---- Count interactions
# input <- list.files("bamFiles/", pattern = "_trimmed.h5", full.names = TRUE)
# names(input) <- str_remove_all(basename(input), "_R1_2.hicup_trimmed.h5")

# mergePairs(files=input, "merged.h5")

# The promoters are from the cHiC paper suppl. data. Promoters with two baits
# were merged as promoters with two baits had twice the number of interactions
# (one for each bait)
# promoters <- import("promoters_hg38_fixed.bed")
# promoters <- unstrand(promoters)
# mcols(promoters) <- NULL
# 
# enhancers <- import("../ChIP/enhancers.bed")
# 
# enhancer <- reduce(resize(enhancers, fix="center", width=1e4))
# viewpoints <- connectCounts(files = input, param = hs.param, regions = promoters, second.regions = c(enhancer))
# 
# 
# ## ---- Filter out extra regions
# filterFn <- function(vps)
# {
#   ave.ab <- aveLogCPM(asDGEList(vps))
#   trended <- filterTrended(vps)
#   
#   count.keep <- ave.ab >= aveLogCPM(5, lib.size=mean(vps$totals))
#   trend.keep <- trended$abundances > trended$threshold + log2(2)
#   final.keep <- count.keep & trend.keep
#   
#   vps[final.keep]
# }
# viewpoints <- filterFn(viewpoints)
# 
# save(viewpoints, file = "selectedInteractions.RData")
# 
# 
# ### Figure of genomes
# chroms <- as(seqinfo(BSgenome.Hsapiens.UCSC.hg38), "GRanges") %>% keepStandardChromosomes(pruning.mode = "coarse")
# 
# hs.frag <- cutGenome(BSgenome.Hsapiens.UCSC.hg38, "AAGCTT", 4)
# hs.param <- pairParam(hs.frag)
# 
# oldPen <- options("scipen")
# options(scipen = 100)
# 
# plotFn <- function(regions, nBins = 1000)
# {
#   op <- par(mar = c(7,7,4,2) + 0.1, las = 2, mgp=c(5,1,0))
#   size <- width(regions)
#   xlab <- sub("chr", "Chromosome ", seqnames(regions))
#   rotPlaid(file = "merged.h5", param = hs.param, region = regions, width = size/nBins, col = "red", max.count = 20, xlab = xlab, ylab = "Gap", axes = FALSE)
#   breaks <- pretty(c(1, size))
#   axis(side = 1, at = breaks, labels = paste(breaks/10^6, "Mbp"))
#   axis(side = 2, at = breaks, labels = paste(breaks/10^6, "Mbp"))
#   par(op)
# }
# 
# for (i in seq_along(chroms)){
#   png(paste0("figures/", seqnames(chroms[i]), ".png"), width = 15, height = 10, units = "cm", res = 300)
#   plotFn(chroms[i])
#   dev.off()
# }
# options(scipen = oldPen)
# 
# # Wed Aug 29 11:13:55 2018 ------------------------------
# # Figure showing interactions relative to a specific region of interest
# ftoRegion <- GRanges("chr16", IRanges(53764083, 53802889))
# startSite <- start(ftoRegion) - 500000
# endSite <- end(ftoRegion) + 1300000
# bins <- seq(from = startSite, to = endSite, by = 5000)
# surroundingRegion <- GRanges("chr16", IRanges(bins[-length(bins)], bins[-1]))
# 
# hs.frag <- cutGenome(BSgenome.Hsapiens.UCSC.hg38, "AAGCTT", 4)
# hs.param <- pairParam(hs.frag)
# 
# countsOfInterest <- connectCounts(files = "merged.h5", 
#                                   param = hs.param, 
#                                   regions = ftoRegion, 
#                                   second.regions = surroundingRegion
#                                   )
# 
# selector <- function(x) x[x$is.second]
# secondRegion <- c(anchors(countsOfInterest, "second") %>% selector,
#                   anchors(countsOfInterest, "first") %>% selector)
# 
# pD <- data.table(counts = assay(countsOfInterest, "counts") %>% as.vector)
# pD[, position:=(start(secondRegion) + end(secondRegion))/2]
# 
# p1 <- ggplot(pD, aes(x = position, y = counts)) + geom_line(size = 3) +
#   #geom_vline(xintercept = start(ftoRegion), lty = 2) +
#   #geom_vline(xintercept = end(ftoRegion), lty = 2) +
#   theme_cbmr() +
#   theme_void() +
#   ylab("Supporting reads") +
#   scale_x_continuous(name = "chr16", expand = c(0,0))
# 
# size = 1000
# tiff(width=size, height=size/3, units="px", filename="figures/ftoPlot.tiff", compression = "lzw+p", type = "cairo")
# p1
# dev.off()