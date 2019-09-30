options(java.parameters = "-Xmx4096m")  ## memory set to 4 GB
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(diffHic)
require(BSgenome.Hsapiens.UCSC.hg38)
require(edgeR)
library(csaw)
library(magrittr)
library(data.table)
library(openxlsx)
library(stringr)
library(GenomicRanges)
load("../selectedInteractions.RData")

runEdgeR <- function(data)
{
  colnames(data) <- c("P6_Ctrl","P6_Palm","P6_TNFa","P8_Ctrl","P8_Palm","P8_TNFa")
  colData(data)$Passage <- str_split(colnames(data), "_", simplify = TRUE)[, 1]
  colData(data)$Treatment <- str_split(colnames(data), "_", simplify = TRUE)[, 2]

  ab <- aveLogCPM(asDGEList(data))
  o <- order(ab)
  adj.counts <- cpm(asDGEList(data), log=TRUE)
  mval <- rowMeans(adj.counts[,c(1:3)])-rowMeans(adj.counts[,c(4:6)])
  smoothScatter(ab, mval, xlab="A", ylab="M")
  fit <- loessFit(x=ab, y=mval)
  lines(ab[o], fit$fitted[o], col="red")

  normfacs <- normOffsets(data)
  normfacs2 <- normOffsets(data, type="loess", se.out=TRUE)
  
  design <- model.matrix(~ 0 + colData(data)$Passage + colData(data)$Treatment)
  colnames(design) <- c("P6", "P8", "Palmitate", "TNFa")
  
  y <- asDGEList(data, group = colData(data)$Treatment)
  y$offset <- assay(normfacs2, "offset")
  
  ind <- filterByExpr(y)
  y <- y[ind, ]
  data <- data[ind, ]

  y <- estimateDisp(y, design = design)
  plotBCV(y)
  
  fit <- glmQLFit(y, design, robust=TRUE)
  plotQLDisp(fit)
  
  result <- list(
    Palmitate = glmQLFTest(fit, coef=3),
    TNFa = glmQLFTest(fit, coef=4),
    Difference = glmQLFTest(fit, contrast = c(0, 0, 1, -1)),
    Average = glmQLFTest(fit, contrast = c(0, 0, 0.5, 0.5)),
    Passage = glmQLFTest(fit, contrast = c(-1, 1, 0, 0))
  )
  clustered <- clusterPairs(data, tol=1, upper=1e6)
  clusteredResult <- lapply(result, . %>% extract2("table") %>% 
                              combineTests(ids = clustered$indices[[1]]) %>% 
                              as.data.table() %>%
                              setkey("PValue")
                            )
  list(
    allResults = result,
    clusteredResults = clusteredResult,
    viewpoints = data
  )
}

edgeRres <- runEdgeR(viewpoints)
save(edgeRres, file = "edgeRresults.RData")
# Make a list of what interacts with what
# First we need to translate to ensembl IDs
promoters <- import("../promoters_hg38_fixed.bed.gz") %>%
  as.data.table
promoters[, name:=str_remove(name, "[+-]_")]
promoters[, name:=str_split(name, pattern = ",")]

ensembl <- useMart('ENSEMBL_MART_ENSEMBL', 
                   host = "feb2014.archive.ensembl.org", 
                   dataset = 'hsapiens_gene_ensembl') %>%
  getBM(attributes= c("ensembl_gene_id",
                      "entrezgene",
                      "external_transcript_id",
                      "hgnc_symbol",
                      "description"
  ), 
  mart = .) %>% 
  data.table
setkey(ensembl, "external_transcript_id")
promoters[, ensIDs:=apply(promoters, 1, function(x){ensembl[x$name, unique(ensembl_gene_id), nomatch = FALSE]})]
promoters[unlist(lapply(promoters$ensIDs, length)) != 0, ]

# Remove all entries in promoter without an ENSEMBL ID
promoters <- promoters[unlist(lapply(promoters$ensIDs, length)) != 0, ]
promoters[, name:=NULL]
promoters[, score:=NULL]
promoters <- as(promoters, "GRanges")

enhancers <- import("../../1.ChIP/enhancers.bed")

promoters2enhancers <- function(hicRes, prom, enh)
{
  vp <- hicRes$viewpoints
  fst <- anchors(vp, type="first")
  fst$hicID <- seq_along(fst)
  snd <- anchors(vp, type="second")
  snd$hicID <- seq_along(snd)
  
  addToGR <- function(anchor, features, ids)
  {
    ind <- findOverlaps(anchor, features)
    out <- tapply(ids[subjectHits(ind)], queryHits(ind), unlist)
    list(inds = unique(queryHits(ind)), ids = out)
  }
  fst$promoterIDs <- ""
  ind <- addToGR(fst, prom, prom$ensIDs)
  fst$promoterIDs[ind$inds] <- ind$ids
  
  snd$promoterIDs <- ""
  ind <- addToGR(snd, prom, prom$ensIDs)
  snd$promoterIDs[ind$inds] <- ind$ids
  
  fst$enhancerIDs <- ""
  ind <- addToGR(fst, enh, enh$name)
  fst$enhancerIDs[ind$inds] <- ind$ids
  
  snd$enhancerIDs <- ""
  ind <- addToGR(snd, enh, enh$name)
  snd$enhancerIDs[ind$inds] <- ind$ids
  
  list(first = fst, second = snd)
}

promoterEnhancerAnchors <- promoters2enhancers(edgeRres, promoters, enhancers)
save(promoterEnhancerAnchors, file = "annotatedAnchors.RData")

unlistAndMerge <- function(gr1, gr2, byStr, featureStr)
{
  gr1 <- gr1 %>% as.data.table()
  gr1 <- gr1[, .(promoterIDs = unlist(promoterIDs)), by = "hicID"]
  gr1 <- gr1[promoterIDs != "", ]
  
  gr2 <- gr2 %>% as.data.table()
  gr2 <- gr2[, .(enhancerIDs = unlist(enhancerIDs)), by = "hicID"]
  gr2 <- gr2[enhancerIDs != "", ]
  
  out <- merge(gr1, gr2, by = "hicID")
  out <- out[, list(hicIDs = list(unique(hicID)), 
                    feature = list(unique(get(featureStr)))), by = byStr]
  out <- out[get(byStr) != ""]
  out
}

makeMaps <- function(anchors, byStr, featureStr)
{
  fst2snd <- unlistAndMerge(anchors$first, anchors$second, byStr, featureStr)
  snd2fst <- unlistAndMerge(anchors$second, anchors$first, byStr, featureStr)
  
  mapDT <- rbind(fst2snd, snd2fst)
  mapDT[, .(hicIDs = list(unique(unlist(hicIDs))), feature = list(unique(unlist(feature)))), by = byStr]
}

promoters2enhancers <- makeMaps(promoterEnhancerAnchors, byStr = "promoterIDs", featureStr = "enhancerIDs")
setnames(promoters2enhancers, "feature", "enhancerID")
enhancers2promoters <- makeMaps(promoterEnhancerAnchors, byStr = "enhancerIDs", featureStr = "promoterIDs")
setnames(enhancers2promoters, "feature", "promoterID")
save(promoters2enhancers, enhancers2promoters, file = "featureMaps.RData")
