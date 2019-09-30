library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(SummarizedExperiment)
library(reactome.db)
library(GO.db)
library(clusterProfiler)
library(openxlsx)

load("edgerDataObjects.RData")
# Camera enrichment analysis

cameraFn <- function(contrast, index, y, design)
{
  out <- camera(y, index = index, design = design, contrast = contrast) %>%
    data.table(keep.rownames = TRUE) %>% 
    setnames("rn", "GOID") %>%
    extract(!is.na(PValue))
  out[order(PValue)]
}

#### Convert from ENSEMBL to ENTREZID
conv <- bitr(rownames(y$counts), fromType='ENSEMBL', toType='ENTREZID', OrgDb = "org.Hs.eg.db") %>%
  data.table(key = "ENSEMBL")
rownames(y) <- conv[rownames(y), ENTREZID, mult = "first"]


#### Gene Ontology (only BP)
keysGO <- keys(GO.db)
termGO <- AnnotationDbi::select(GO.db, keys=keysGO, columns=c("TERM", "ONTOLOGY")) %>% data.table
termGO <- termGO[ONTOLOGY == "BP"]
termGO[, ONTOLOGY:=NULL]

cyt.go.genes <- as.list(org.Hs.egGO2ALLEGS)
cyt.go.genes <- cyt.go.genes[names(cyt.go.genes) %in% termGO$GOID]
cyt.go.genes <- Filter(. %>% length %>% is_greater_than(4), cyt.go.genes) # Remove small categories
cyt.go.genes <- Filter(. %>% length %>% is_less_than(501), cyt.go.genes) # Remove large categories

cameraGO <- apply(ctrsts, 2, cameraFn, index = cyt.go.genes, y = y, design = design)
cameraGO <- lapply(cameraGO, merge, termGO, by = "GOID", nomatch = FALSE)
cameraGO <- lapply(cameraGO, extract, order(PValue, decreasing = FALSE)) 

#### Reactome
reactomeName <- as.list(reactomePATHID2NAME)
reactomeName <- data.table(GOID = names(reactomeName), TERM = reactomeName)
reactomeList <- as.list(reactomePATHID2EXTID)
reactomeList <- Filter(. %>% length %>% is_greater_than(4), reactomeList) # Remove small categories
reactomeList <- Filter(. %>% length %>% is_less_than(501), reactomeList) # Remove small categories

cameraReactome <- apply(ctrsts, 2, cameraFn, index = reactomeList, y = y, design = design)
cameraReactome <- lapply(cameraReactome, merge, reactomeName, by = "GOID", nomatch = FALSE)
cameraReactome <- lapply(cameraReactome, extract, order(PValue, decreasing = FALSE)) 

write.xlsx(cameraGO, "edgeR_results/cameraGO.xlsx", asTable = TRUE)
write.xlsx(cameraReactome, "edgeR_results/cameraReactome.xlsx", asTable = TRUE)
