library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(SummarizedExperiment)
library(pheatmap)
library(scales)
require(factoextra)
library(RColorBrewer)
library(openxlsx)
require(reshape2)
require(gridExtra)
library(vegan)
require(dplyr)
require(zoo)

exprData <- read.xlsx("../qcFigures/cpmData.xlsx", sheet = 1) %>%
  data.table()
exprData[, c("SYMBOL", "GENENAME"):=NULL]
ensIDs <- exprData[, 1] %>% unlist
exprData <- as.matrix(exprData[, -1])
rownames(exprData) <- unname(ensIDs)

colDat <- data.frame(tstrsplit(colnames(exprData), "_"))
colnames(colDat) <- c("Passage", "Treatment", "Replicate")
levels(colDat$Treatment) <- c("Control", "Palmitate", "TNFa")
levels(colDat$Replicate) <- c("Rep1", "Rep2")
colDat <- lapply(colDat, as.character) %>% data.frame(stringsAsFactors = FALSE)
sumExp <- SummarizedExperiment(list("cpmVals" = exprData),
                               colData = colDat
                               )

formatForHeatmap <- function(sumExp, selectedGenes = NULL, selectedSamples = NULL)
{
  if (!is.null(selectedGenes)){
    sumExp <- sumExp[selectedGenes, ]
  }
  if (!is.null(selectedSamples)){
    sumExp <- sumExp[, selectedSamples]
  }
  rawLogCPM <- log(assays(sumExp)[["cpmVals"]] + 0.25)
  
  scaled <- t(scale(t(rawLogCPM), center = TRUE, scale = TRUE))
  # Some genes do not vary across samples (usually 0 counts). Remove them:
  hasNA <- apply(scaled, 1, anyNA)
  scaled <- scaled[!hasNA, ]
  
  colAnnot <- colData(sumExp)
  
  list(scaled = scaled,
       colAnnot = colAnnot
  )
}

clusterTester <- function(heatmapData, testName)
{
  distMat <- as.dist(1 - cor(t(heatmapData$scaled))/2)
  # Calculate optimal k
  pdf(paste0("./heatmaps/", testName, "_optimal_k.pdf"))
  d    <- heatmapData$scaled
  kmax <- 15
  nbc  <- list()
  nbcMethod <- c("silhouette", "wss", "gap_stat")
  for (i in 1:length(nbcMethod)){
    nbc[[i]] <- fviz_nbclust(d, kmeans, diss = distMat, method = nbcMethod[i], k.max = kmax)
    print(nbc[[i]])
  }
  dev.off()
}

edgerResultsFile <- "../edgeR_results/RNA.xlsx"
sheets <- getSheetNames(edgerResultsFile)
names(sheets) <- c("Palmitate", "TNFa", "Replicate")

signifGenes <- lapply(sheets, . %>% 
                        read.xlsx(xlsxFile = edgerResultsFile) %>% 
                        data.table %>% 
                        extract(FDR < 0.05, ENSEMBL))

hmData <- list(
  all = formatForHeatmap(sumExp),
  Palmitate = formatForHeatmap(sumExp, selectedGenes = signifGenes$Palmitate),
  TNFa = formatForHeatmap(sumExp, selectedGenes = signifGenes$TNFa),
  Replicate = formatForHeatmap(sumExp, selectedGenes = signifGenes$Replicate)
)

# clusterTester(hmData$Palmitate, "Palmitate")
# clusterTester(hmData$TNFa, "TNFa")
# clusterTester(hmData$Replicate, "Replicate")


# hmCols <- cet_pal("cyclic_mygbm_30-95_c78_n256", n = 256)
# hmCols[floor(c(6, 10, 18, 22)/24 * 256)]

makeColours <- function(heatmapData)
{
  baseCols <- brewer.pal(7, "RdYlBu") %>% rev
  dataMat <- heatmapData$scaled
  quantiles <- quantile(dataMat, probs = c(0.10, 0.90))
  highestVal <- max(abs(quantiles))
  hmBreaks <- c(min(c(min(dataMat), -highestVal)) - 0.1, seq(-highestVal, highestVal, length.out = 256), max(c(max(dataMat), highestVal)) + 0.1)
  hmCols <- gradient_n_pal(baseCols)(seq(0, 1, length=length(hmBreaks)-2))
  hmCols <- c(baseCols[1], hmCols, baseCols[7])
  
  annotColour <- list(
    Passage = c("P5" = "#66c2a4", 
                "P6" = "#2ca25f"
                  ),
    Treatment = c(Control = "#1b9e77",
                  Palmitate = "#d95f02",
                  TNFa = "#7570b3"),
    Replicate = c(Rep1 = "#a6cee3",
                  Rep2 = "#1f78b4")
  )
  annotColour$Treatment <- annotColour$Treatment[names(annotColour$Treatment) %in% heatmapData$colAnnot$Treatment]
  list(hmCols = hmCols, annotColour = annotColour)
}

heatmapfn <- function(heatmapData, name)
{
  colours <- makeColours(heatmapData)
  
  anot <- as.data.table(heatmapData$colAnnot)
  anot[, n:=1:.N]
  anot[, Treatment:=factor(Treatment, levels = c("Control", "Palmitate", "TNFa"))]
  anot[, Passage:=factor(Passage, levels = c("P5", "P6", "P8"))]
  anot[, Replicate:=factor(Replicate, levels = c("Rep1", "Rep2"))]
  setkey(anot, Treatment, Passage, Replicate)
  anot[, n2:=1:.N]
  setkey(anot, n)
  newOrder <- anot$n2
  
  clustering_callback = function(hclObj, dat)
  {
    if (length(hclObj$order) == nrow(heatmapData$colAnnot)){
      hclObj <- reorder(hclObj, newOrder, agglo.FUN = "mean")
    }
    hclObj
  }
  
  pheatmap(heatmapData$scaled, 
           scale = "none", 
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           clustering_method = 'ward.D2',
           annotation_col = data.frame(heatmapData$colAnnot, stringsAsFactors = FALSE),
           annotation_names_row = FALSE,
           annotation_names_col = FALSE, 
           cutree_rows = 2,
           annotation_colors = colours$annotColour,
           clustering_callback = clustering_callback,
           color = colours$hmCols, 
           show_rownames = FALSE,
           show_colnames = FALSE,
           drop_levels = TRUE,
           filename = paste0("heatmaps/", name, "_pheatmap.png")
  ) 
}

heatmapfn(hmData$Replicate, "Replicate")
heatmapfn(hmData$Palmitate, "Palmitate")
heatmapfn(hmData$TNFa, "TNFa")
heatmapfn(hmData$all, "all")



hmData_subsets <- list(
  all_palmitate = formatForHeatmap(sumExp, 
                                   selectedSamples = colnames(sumExp)[sumExp$Treatment %in% c("Control", "Palmitate")]),
  all_TNFa = formatForHeatmap(sumExp, 
                              selectedSamples = colnames(sumExp)[sumExp$Treatment %in% c("Control", "TNFa")]),
  Palmitate = formatForHeatmap(sumExp, 
                               selectedGenes = signifGenes$Palmitate, 
                               selectedSamples = colnames(sumExp)[sumExp$Treatment %in% c("Control", "Palmitate")]),
  TNFa = formatForHeatmap(sumExp, 
                          selectedGenes = signifGenes$TNFa, 
                          selectedSamples = colnames(sumExp)[sumExp$Treatment %in% c("Control", "TNFa")])
)
heatmapfn(hmData_subsets$all_palmitate, "all_subset_Palmitate")
heatmapfn(hmData_subsets$all_TNFa, "all_subset_TNFa")
heatmapfn(hmData_subsets$Palmitate, "Palmitate_subset")
heatmapfn(hmData_subsets$TNFa, "TNFa_subset")
