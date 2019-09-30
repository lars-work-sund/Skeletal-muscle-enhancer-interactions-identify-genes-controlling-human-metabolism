library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(SummarizedExperiment)
library(gridExtra)
library(gtable)
library(grid)

xlsxFile <- "../edgeR_results/RNA.xlsx"
getSheetNames(xlsxFile)
logFCdata <- list(
  Palmitate = read.xlsx(xlsxFile, "Palmitate") %>% 
    data.table() %>% 
    setnames(c("logFC", "PValue", "FDR"), 
             c("logFC_Palmitate", "PValue_Palmitate", "FDR_Palmitate")),
  TNFa = read.xlsx(xlsxFile, "TNFa") %>% 
    data.table() %>%
    setnames(c("logFC", "PValue", "FDR"), 
             c("logFC_TNFa", "PValue_TNFa", "FDR_TNFa"))
)
pD <- merge(logFCdata$Palmitate[, .(ENSEMBL, logFC_Palmitate, PValue_Palmitate, FDR_Palmitate)],
            logFCdata$TNFa[, .(ENSEMBL, logFC_TNFa, PValue_TNFa, FDR_TNFa)])

pD[, cor.test(logFC_Palmitate, logFC_TNFa, method = "spearman")]


pBase <- ggplot(pD[order(PValue_Palmitate, decreasing = TRUE), ], 
                aes(logFC_Palmitate, logFC_TNFa)) + 
  geom_point(size = 2, colour = "grey50", alpha = 0.5) +
  geom_smooth(method='lm',formula=y~x, se = FALSE, colour = "black", linetype = 2) +
  theme(legend.justification = c(0, 0.5)) +
  theme_bw(base_size = 30) +
  scale_x_continuous(name = expression(log[2](FC[paste("Palmitate")]))) +
  scale_y_continuous(name = expression(log[2](FC[paste("TNF", alpha)]))) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_minimal(base_size = 24)


size = 1000
tiff(width=size, height=size, units="px", filename="logFCs/RNA_logFC.tiff", compression = "lzw+p", type = "cairo")
pBase %>%
  print
dev.off()

