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

load("edgeRresults.RData")
load("featureMaps.RData")

scorePromoters <- function(tableStats, promoterMap)
{
  tableStats <- data.table(tableStats)
  scoreFn <- function(ind) tableStats[unlist(ind), sum(-log10(PValue) * logFC)]
  out <- promoterMap[, .(score = scoreFn(hicIDs)), by = "promoterIDs"]
  out[order(score, decreasing = TRUE), ]
}

promoterScores <- list(
  palmitate = scorePromoters(edgeRres$allResults$Palmitate$table, promoters2enhancers),
  tnfa = scorePromoters(edgeRres$allResults$TNFa$table, promoters2enhancers),
  passage = scorePromoters(edgeRres$allResults$Passage$table, promoters2enhancers)
)

# Camera might be a bad fit for this analysis, it assumes that correlations happens within gene sets, 
# here correlations happens withing genes close to each oter.

ensembl <- useMart('ENSEMBL_MART_ENSEMBL', 
                   host = "feb2014.archive.ensembl.org", 
                   dataset = 'hsapiens_gene_ensembl') %>%
  getBM(attributes= c("ensembl_gene_id",
                      "hgnc_symbol",
                      "entrezgene",
                      "description"
  ), 
  mart = .) %>% 
  data.table
setkey(ensembl, ensembl_gene_id)

promoterScores <- lapply(promoterScores, . %>% ensembl[., mult = "first"])
write.xlsx(promoterScores, file = "promoterScores.xlsx", asTable = TRUE)
