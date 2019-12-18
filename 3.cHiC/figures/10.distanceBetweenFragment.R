library(diffHic)
library(ggplot2)
load("../selectedInteractions.RData")
data <- viewpoints

labFn <- function(x)
{
  out <- x
  x <- na.omit(x)
  suffix <- character(length(x))
  suffix[abs(x) < 10^3] <- "bp"
  suffix[abs(x) >= 10^3] <- "kb"
  suffix[abs(x) >= 10^6] <- "Mb"
  suffix[abs(x) >= 10^9] <- "Gb"
  
  ind1 <- abs(x) >= 10^3 & abs(x) < 10^6
  ind2 <- abs(x) >= 10^6 & abs(x) < 10^9
  ind3 <- abs(x) >= 10^9
  x[ind1] <- round(x[ind1]/10^3, 2)
  x[ind2] <- round(x[ind2]/10^6, 2)
  x[ind3] <- round(x[ind3]/10^9, 2)
  out[!is.na(out)] <- paste(x, suffix)
  out
}

dist <- distance(anchors(data)[[2]], anchors(data)[[1]])
dens <- density(na.omit(log10(dist)))

modeDist <- 10^(dens$x[which.max(dens$y)])
medianDist <- median(dist, na.rm = TRUE)
meanDist <- mean(dist, na.rm = TRUE)

ggplot(data.table(dist), aes(x = dist)) + 
  geom_histogram(bins = 45) +
  scale_x_continuous(labels = labFn, trans = "log10", name = "Distance between fragments") +
  scale_y_continuous(name = "Count") +
  annotation_logticks(sides = "b") +
  theme_cbmr() +
  theme(plot.margin=unit(c(5,10,5,5),"mm")) +
  geom_vline(xintercept = medianDist, lty = 3)

ggsave(filename = "distanceBetweenFragments.pdf")
meanDist
