library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(cbmR)
library(SummarizedExperiment)

promoters <- import("promoters_hg38.bed.gz")
promoters$nameUniq <- sub("^(\\+|\\-)_", "", promoters$name)

# Remove the ultra conserved elements (uce)
uceS <- promoters[promoters$nameUniq == "uce"]
promoters <- promoters[promoters$nameUniq != "uce"]

# Seperate into genes tageted once and twice
dupGenes <- promoters$nameUniq[duplicated(promoters$nameUniq)]
singleGenes <- setdiff(promoters$nameUniq, dupGenes)

promotersSingle <- promoters[promoters$nameUniq %in% singleGenes]
promotersDouble <- promoters[promoters$nameUniq %in% dupGenes]

# Check that every gene appears exactly twice
promotersDouble$nameUniq %>% table %>% equals(2) %>% all

### Merge the two ends
promotersDouble %<>% 
  as.data.table %>%
  extract(, .(seqnames = seqnames[[1]],
              start = min(start),
              end = max(end),
              name = name[[1]],
              score = 0),
          by = "nameUniq") %T>%
  setcolorder(c(2:6,1)) %>%
  as("GRanges")

promotersFixed <- c(uceS, promotersSingle, promotersDouble)
export(promotersFixed, con = "../promoters_hg38_fixed.bed.gz")