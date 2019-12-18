library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(cbmR)
library(clipr)
library(openxlsx)
library(clusterProfiler)
library(GO.db)
library(org.Hs.eg.db)
library(viridisLite)
library(cowplot)
library(export)
extract <- magrittr::extract
goToFile()

load("../edgerDataObjects.RData")
cyt.go.genes <- as.list(org.Hs.egGO2ALLEGS)



addFractionSignificant <- function(type, edgerRes, cyt.go.genes)
{
  edgerRes <- copy(edgerRes[[type]])
  conv <- bitr(edgerRes$ENSEMBL, fromType='ENSEMBL', toType='ENTREZID', OrgDb = "org.Hs.eg.db") %>%
    data.table(key = "ENSEMBL")
  edgerRes[, ENTREZ:=conv[ENSEMBL, ENTREZID, mult = "first"]]
  setkey(edgerRes, ENTREZ)
  
  goRes <- read.xlsx("edgeR_results/cameraGO.xlsx", type) %>% setDT
  down <- goRes[Direction == "Down"][1:10]
  up <- goRes[Direction == "Up"][1:10]
  
  fracSignif <- function(x){
    fdrCutoff <- 0.05
    edgerRes[cyt.go.genes[x], sum(FDR < fdrCutoff)/.N, nomatch = FALSE]
  }
  up[, fractionSignificant:=lapply(GOID, fracSignif) %>% unlist]
  down[, fractionSignificant:=lapply(GOID, fracSignif) %>% unlist]
  rbind(up, down)
}

sheets <- c("Palmitate", "TNFa")
pDs <- lapply(sheets, addFractionSignificant, edgerRes = annotated, cyt.go.genes)
names(pDs) <- sheets

minusLog10 <- scales::trans_new("minusLog10",
                                transform = function(x){-log10(x)},
                                inverse = function(x){10^(-x)},
                                breaks = function(x){10^(-scales::pretty_breaks()(-log10(x)))},
                                format = function(x){scales::math_format()(log10(x))},
                                domain = c(0, 1))
limits <- do.call("rbind", pDs)[, .(PValue = range(PValue), NGenes =range(NGenes), fractionSignificant = range(fractionSignificant))]

colScale <- scale_colour_gradient(name = "P-value", trans = minusLog10, low = "#9ebcda", high = "#4d004b", limits = rev(limits$PValue))
sizeScale <- scale_size_continuous(name = "Genes in\n term", trans = "log10", limits = limits$NGenes)
xScale <- scale_x_continuous("Percent of genes significant", labels = scales::percent_format(accuracy = 1), limits = limits$fractionSignificant) 

pDs$Palmitate[, Comparison:="Palmitate"]
pDs$TNFa[, Comparison:="TNFa"]

pD <- do.call("rbind", pDs)
pD[, TERM:=ordered(TERM, levels = rev(TERM[order(PValue)]))]
pD[, Direction:=ordered(Direction, levels = c("Up", "Down"))]
bigPlot <- ggplot(pD, aes(x = fractionSignificant, y = TERM, size = NGenes, colour = PValue)) +
  geom_point() +
  colScale +
  sizeScale +
  xScale +
  scale_y_discrete(name = NULL, ) +
  theme_cowplot() +
  facet_wrap( ~ Direction + Comparison, scales = "free", ncol = 2)
ggsave("goFigures/goPlot.pdf", bigPlot, width = 50, height = 20, units = "cm")
inchToCm <- 0.393701
graph2office(bigPlot, file = "goFigures/goPlot.pptx", type = "PPT", width = 50 * inchToCm, height = 20 * inchToCm)

plotFn <- function(pD){
  subFn <- function(direction){
    ggplot(pD[Direction == direction], aes(x = fractionSignificant, y = TERM, size = NGenes, colour = PValue)) +
      geom_point() +
      colScale +
      sizeScale +
      xScale +
      scale_y_discrete(name = NULL, limits = pD[Direction == direction, rev(TERM)], drop = TRUE) +
      theme_cowplot()
  }
  list(upPlot = subFn("Up"),
       downPlot = subFn("Down"))
}

tnfaPlots <- plotFn(pDs$TNFa)
palmPlots <- plotFn(pDs$Palmitate)

ggsave("goFigures/tnfa_goterms_up.pdf", tnfaPlots$upPlot, width = 20, height = 10, units = "cm")
ggsave("goFigures/tnfa_goterms_down.pdf", tnfaPlots$downPlot, width = 20, height = 10, units = "cm")
ggsave("goFigures/palm_goterms_up.pdf", palmPlots$upPlot, width = 20, height = 10, units = "cm")
ggsave("goFigures/palm_goterms_down.pdf", palmPlots$downPlot, width = 20, height = 10, units = "cm")
