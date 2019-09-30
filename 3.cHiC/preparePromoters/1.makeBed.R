library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(SummarizedExperiment)
library(rtracklayer)
library(stringr)
library(data.table)
promoters <- import("promoters.fasta.gz")
dt <- data.table(as.data.frame(str_split_fixed(names(promoters), ":|-| ", 4), stringsAsFactors = FALSE))
setnames(dt, c("seqnames", "start", "end", "name"))
dt[, c("strand", "name"):=as.data.frame(str_split_fixed(name, " ", 2), stringsAsFactors = FALSE)]
dt[, c("start", "end"):=list(as.integer(start), as.integer(end))]
dt[, seqnames:=sub("C", "c", seqnames)]

fn <- function(seqnames, start, end)
{
  uniqSeq <- unique(seqnames)
  if (length(uniqSeq) != 1) data.table(seqnames = seqnames, start = start, end = end)
  minStart <- min(start)
  maxEnd <- max(end)
  if (maxEnd - minStart <= 10000)
  {
    list(seqnames = uniqSeq, start = minStart, end = maxEnd)
  } else
  {
    data.table(seqnames = seqnames, start = start, end = end)
  }
}

dt2 <- dt[, fn(seqnames, start, end), by = "name"]

## Regions are rounded to the nearest restriction site. Overlaps should be fine, use the entire set.
gr <- GRanges(dt)
gr$name <- paste(strand(gr), gr$name, sep = "_")

export(gr, con = "promoters.bed.gz")
