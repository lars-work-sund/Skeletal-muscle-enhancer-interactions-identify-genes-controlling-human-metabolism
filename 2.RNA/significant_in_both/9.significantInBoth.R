library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(cbmR)
library(clipr)
library(clusterProfiler)
extract <- magrittr::extract

load("../edgerDataObjects.RData")
signifInBoth <- intersect(annotated$Palmitate[FDR < 0.05, SYMBOL], 
                          annotated$TNFa[FDR < 0.05, SYMBOL])

signifInEither <- union(annotated$Palmitate[FDR < 0.05, SYMBOL], 
                          annotated$TNFa[FDR < 0.05, SYMBOL])

# The same genes are tested in both conditions
allGenes <- annotated$Palmitate[, SYMBOL]


enrichments <- list(
  vsAllGenes = enrichGO(signifInBoth, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", universe = allGenes) %>% as.data.table,
  vsAllSignificant = enrichGO(signifInBoth, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", universe = signifInEither) %>% as.data.table
)

write.xlsx(enrichments, file = "significantInBoth.xlsx")
