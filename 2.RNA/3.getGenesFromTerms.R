library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(clusterProfiler)
library(SummarizedExperiment)
library(openxlsx)

# Mon Jun 11 12:34:03 2018 ------------------------------
# Extract genes related to interesting GO terms
load("edgerDataObjects.RData")

#### Gene Ontology (only BP)
keysGO <- keys(GO.db)
termGO <- AnnotationDbi::select(GO.db, keys=keysGO, columns=c("TERM", "ONTOLOGY")) %>% data.table
termGO <- termGO[ONTOLOGY == "BP"]
termGO[, ONTOLOGY:=NULL]

cyt.go.genes <- as.list(org.Hs.egGO2ALLEGS)
cyt.go.genes <- cyt.go.genes[names(cyt.go.genes) %in% termGO$GOID]
cyt.go.genes <- Filter(. %>% length %>% is_greater_than(4), cyt.go.genes) # Remove small categories
cyt.go.genes <- Filter(. %>% length %>% is_less_than(501), cyt.go.genes) # Remove large categories

getGenesFromTerm <- function(term, genes)
{
  conv <- bitr(genes$ENSEMBL, fromType='ENSEMBL', toType='ENTREZID', OrgDb = "org.Hs.eg.db") %>%
    data.table(key = "ENSEMBL")
  genes[, ENTREZID:=conv[ENSEMBL, ENTREZID, mult = "first"]]
  termOfInterest <- termGO[TERM == term]
  genesOfInterest <- cyt.go.genes[[termOfInterest$GOID]]
  if (is.null(genesOfInterest)) return(NA)
  out <- genes[genesOfInterest, on = "ENTREZID", nomatch = FALSE]
  out[, ENTREZID:=NULL]
  out <- unique(out)
  out <- out[order(PValue, decreasing = FALSE), ]
  out
}
palmTerms <- c("positive regulation of inflammatory response", "regulation of lipid metabolic process", "muscle filament sliding", "actin-myosin filament sliding", "nucleosome assembly", "regulation of fatty acid oxidation", "regulation of lipid biosynthetic process")
names(palmTerms) <- palmTerms
TNFaTerms <- c("positive regulation of inflammatory response", "muscle filament sliding", "actin-myosin filament sliding", "regulation of insulin-like growth factor receptor signaling pathway", "type B pancreatic cell proliferation", "regulation of lipid biosynthetic process")
names(TNFaTerms) <- TNFaTerms

palmGenes <- lapply(palmTerms, getGenesFromTerm, genes = copy(annotated$Palm))
TNFaGenes <- lapply(TNFaTerms, getGenesFromTerm, genes = copy(annotated$TNFa))

write.xlsx(palmGenes, file = "edgeR_results/palmTermsGenes.xlsx", asTable = TRUE)
write.xlsx(TNFaGenes, file = "edgeR_results/TNFaTermsGenes.xlsx", asTable = TRUE)



terms <- c("acute inflammatory response",
           "muscle filament sliding",
           "nucleosome assembly",
           "protein targeting to ER",
           "regulation of lipid metabolic process",
           "type I interferon signaling pathway",
           "cellular response to tumor necrosis factor")
names(terms) <- terms

palmGenes <- lapply(terms, getGenesFromTerm, genes = copy(annotated$Palm))
TNFaGenes <- lapply(terms, getGenesFromTerm, genes = copy(annotated$TNFa))

write.xlsx(palmGenes, file = "edgeR_results/palmTermsGenes2.xlsx", asTable = TRUE)
write.xlsx(TNFaGenes, file = "edgeR_results/TNFaTermsGenes2.xlsx", asTable = TRUE)
