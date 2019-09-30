library(data.table)
library(magrittr)
library(edgeR)
library(stringr)
library(clusterProfiler)
library(openxlsx)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library(org.Hs.eg.db)
library(SummarizedExperiment)

makeSumExp <- function(file)
{
  counts <- fread(file)
  countsMat <- counts[, -1] %>% as.matrix
  rownames(countsMat) <- counts$Geneid
  
  countsMat <- countsMat[!rowSums(countsMat) == 0, ]
  expSetup <- data.table(id = colnames(countsMat))
  expSetup[, c('Passage', 'Treatment', 'Replicate'):=tstrsplit(id, "_")]
  expSetup[, Passage:=factor(Passage, levels = c("P5", "P6"))]
  expSetup[, Treatment:=factor(Treatment, levels = c("Ctrl", "Palm", "TNFa"))]
  expSetup[, Replicate:=factor(Replicate, levels = c("Rep1", "Rep2"))]
  
  cDat <- as.data.frame(expSetup)
  rownames(cDat) <- cDat$id
  cDat$id <- NULL
  
  SummarizedExperiment(assays = list(counts = countsMat), colData = cDat)
}


colScale <- scale_colour_manual(
  values = c(Ctrl = "#1b9e77",
             Palm = "#d95f02",
             TNFa = "#7570b3"),
  labels = c(Ctrl = "Control",
             Palm = "Palmitate",
             TNFa = "TNFa")
)

lonzaEnhancerRNA <- makeSumExp("featureCounts.csv.gz")

# Prepare the data: Filter and perform QC plotting
source("../supportFunctions.R", local = TRUE)

y <- rnaPrepare(sumExp = lonzaEnhancerRNA,
                folder = "qcFigures", 
                group = "Treatment",
                mdsLabel = "Passage", 
                mdsColour = "Treatment", 
                colScale = colScale)

# Next create the design matrix
y$samples$subgroup <- with(y$samples, paste(Passage, Treatment, sep = "_"))
design <- model.matrix(~ 0 + subgroup + Replicate, data = y$samples)
colnames(design) <- str_replace_all(colnames(design), "Treatment|Passage|subgroup", "")
colnames(design) <- str_replace(colnames(design), "Replicate", "")
colnames(design) <- str_replace(colnames(design), ":", "_")
# Fit the model
y <- estimateDisp(y, design = design, robust = TRUE) %T>%
  redirectToFile(file.path("qcFigures", "BCV.pdf"), plotBCV)

efit <- glmQLFit(y, design = design, robust = TRUE) %T>%
  redirectToFile(file.path("qcFigures", "QLDispersion.pdf"), plotQLDisp)

colnames(design)
ctrsts <- makeContrasts(
  Palmitate = (P5_Palm + P6_Palm)/2 - (P5_Ctrl + P6_Ctrl)/2,
  TNFa = (P5_TNFa + P6_TNFa)/2 - (P5_Ctrl + P6_Ctrl)/2,
  Replicate = Rep2,
  levels = design
)

dgeResults <- apply(ctrsts, 2, . %>% 
                      glmQLFTest(glmfit = efit, contrast = .) %>% 
                      topTags(n = Inf, p.value = 1) %>% 
                      extract2("table") %>% 
                      as.data.table(keep.rownames = TRUE))

annotated <- lapply(dgeResults, annotateFromEns)

dir.create("edgeR_results", showWarnings = FALSE)
write.xlsx(annotated, file = "edgeR_results/RNA.xlsx", asTable = TRUE)

save(y, efit, design, ctrsts, annotated, file = "edgerDataObjects.RData")
