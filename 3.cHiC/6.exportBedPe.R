library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(SummarizedExperiment)
library(clipr)
library(GenomicInteractions)

load("secondaryStatistics/edgeRresults.RData")
viewpoints <- edgeRres$viewpoints
interactionObject <- GInteractions(anchor1 = anchors(viewpoints)$first,
              anchor2 = anchors(viewpoints)$second)

scoreVector <- cpmByGroup(assay(viewpoints), log = TRUE, group = rep("a", 6))
interactionObject$scoreCol <- as.vector(scoreVector)
mcols(interactionObject)
export.bedpe(GIObject = interactionObject, 
             fn = "selectedInteractions.bedpe",
             score = "scoreCol"
            )
