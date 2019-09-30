library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(SummarizedExperiment)
library(clipr)
library(openxlsx)

readXlsx <- function(xlsxFile)
{
  sheets <- getSheetNames(xlsxFile)
  names(sheets) <- sheets
  lapply(sheets, . %>% read.xlsx(xlsxFile = xlsxFile) %>% data.table)
}

pValues <- readXlsx("p_values_all.xlsx")
pValues <- lapply(pValues, data.table)
for (i in names(pValues)){
  pValues[[i]][, c("Tissue", "Diet"):=tstrsplit(i, "_")]
}
lapply(pValues, setnames, "X1", "Phenotype") %>% invisible
pValues <- do.call("rbind", pValues)
pValues[, Macf1:=as.numeric(Macf1)]
pValues <- melt(pValues, id.vars = c("Diet", "Tissue", "Phenotype"), value.name = "Pvalue", variable.name = "Gene")
pValues[, fdr:=p.adjust(Pvalue, method = "fdr"), by = c("Tissue", "Diet", "Gene")]
out <- dcast(pValues, formula = Phenotype + Tissue + Diet ~ Gene, value.var = "fdr")
write.xlsx(out, "p_values_all_fdr.xlsx")
