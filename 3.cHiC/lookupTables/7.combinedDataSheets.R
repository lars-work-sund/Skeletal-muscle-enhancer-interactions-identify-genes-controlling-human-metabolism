library(openxlsx)
library(biomaRt)
library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(SummarizedExperiment)
library(edgeR)
library(GO.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

load("../secondaryStatistics/edgeRresults.RData")
load("../secondaryStatistics/featureMaps.RData")
load("../secondaryStatistics/annotatedAnchors.RData")

rnaResults <- list(
  Palmitate = read.xlsx("../../2.RNA/edgeR_results/RNA.xlsx", "Palmitate") %>% data.table(key = "ENSEMBL"),
  TNFa = read.xlsx("../../2.RNA/edgeR_results/RNA.xlsx", "TNFa") %>% data.table(key = "ENSEMBL")
)

enhancerResults <- list(
  Palmitate = read.xlsx("../../1.ChIP/enhancerPeaksAnnotated.xlsx", "Palmitate") %>% data.table(key = "rn"),
  TNFa = read.xlsx("../../1.ChIP/enhancerPeaksAnnotated.xlsx", "TNFa") %>% data.table(key = "rn")
)

mergeHiCdata <- function(anchors, edgeRres)
{
  fixAnchor <- function(anchor)
  {
    #anchor$enhancerIDs <- lapply(anchor$enhancerIDs, paste, collapse = ";")
    #anchor$promoterIDs <- lapply(anchor$promoterIDs, paste, collapse = ";")
    as.data.table(anchor)
  }
  
  fst <- fixAnchor(anchors$first)
  setnames(fst, paste0(colnames(fst), "_first"))
  snd <- fixAnchor(anchors$second)
  setnames(snd, paste0(colnames(snd), "_second"))
  
  tab <- edgeRres$table %>% data.table
  setnames(tab, paste0(colnames(tab), "_HiC"))
    
  out <- cbind(fst, snd, edgeRres$table)
  out[, FDR:=p.adjust(PValue)]
  out[, interactionID:=paste0("i", seq_len(nrow(out)))]
  stopifnot(all(out$hicID_first == out$hicID_second))
  out[, c("strand_first", "is.second_first", "original_first", "hicID_first", "strand_second", "nfrags_second", "is.second_second", "hicID_second"):=NULL]
  out
}
mergedHiC <- list(
  Palmitate = mergeHiCdata(promoterEnhancerAnchors, edgeRres$allResults$Palmitate),
  TNFa = mergeHiCdata(promoterEnhancerAnchors, edgeRres$allResults$TNFa),
  Passage = mergeHiCdata(promoterEnhancerAnchors, edgeRres$allResults$Passage)
)

updateHicIDs <- function(x) x[, hicIDs:=lapply(hicIDs, . %>% paste0("i", .))]
updateHicIDs(enhancers2promoters)
setkey(enhancers2promoters, enhancerIDs)
updateHicIDs(promoters2enhancers)
setkey(promoters2enhancers, promoterIDs)
setkey(mergedHiC$Palmitate, interactionID)
setkey(mergedHiC$TNFa, interactionID)
setkey(mergedHiC$Passage, interactionID)
save(enhancers2promoters, promoters2enhancers, mergedHiC, rnaResults, enhancerResults, file = "combinedResults.RData")

numberOfInteractions <- promoters2enhancers[, .(promoterIDs = promoterIDs,
                        nInteractions = lapply(hicIDs, length) %>% unlist, 
                        nEnhancers = lapply(enhancerID, length) %>% unlist)]

conv <- bitr(numberOfInteractions$promoterIDs, fromType = "ENSEMBL", c("SYMBOL", "GENENAME"), "org.Hs.eg.db") %>%
  data.table
numberOfInteractions = merge(numberOfInteractions, conv, by.x = "promoterIDs", by.y = "ENSEMBL", all = TRUE) 
numberOfInteractions <- numberOfInteractions[unique(numberOfInteractions$promoterIDs), on = "promoterIDs", mult = "first"]
write.xlsx(numberOfInteractions, file = "combinedDataSheets/interactionPerPromoter.xlsx")


write.xlsx(enhancers2promoters, file = "combinedDataSheets/enhancers2promoters.xlsx")
write.xlsx(promoters2enhancers, file = "combinedDataSheets/promoters2enhancers.xlsx")

mergeListFn <- function(x)
{
  x[, promoterIDs_first:=lapply(promoterIDs_first, paste, collapse = ", ")]
  x[, enhancerIDs_first:=lapply(enhancerIDs_first, paste, collapse = ", ")]
  x[, promoterIDs_second:=lapply(promoterIDs_second, paste, collapse = ", ")]
  x[, enhancerIDs_second:=lapply(enhancerIDs_second, paste, collapse = ", ")]
  x
}
mergedHiC2 <- lapply(mergedHiC, mergeListFn)
write.xlsx(mergedHiC2, file = "combinedDataSheets/mergedHiC.xlsx")
