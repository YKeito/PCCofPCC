#"/home/yasue/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/DAPSeq2/ATTEEDII_UMAP13PeakIndex_Heatmap.R"
#pacckage----
library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gtools)
#processing data----
before <- proc.time()
#input data----
TAIR10 <- readRDS(file = "/home/yasue/bigdata/yasue/TAIR10_ShortName.rds")
PeakIndex <- readRDS(file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/PeakIndex.rds")
PeakIndex <- PeakIndex %>% select(-contains("PeakCount"))
T.data <- PeakIndex %>% filter(Cluster == "UMAPCluster13")
T.AGI <- data.frame(AGI = permutations(n = length(union(T.data$sAGI, T.data$tAGI)), 
                                       r = 2, v = union(T.data$sAGI, T.data$tAGI), repeats.allowed = T))
T.AGI <- T.AGI %>% mutate(PeakIndex = rep(100, times = nrow(T.AGI)), 
                          sameTFfamily = rep(TRUE, times = nrow(T.AGI)), 
                          PeakMax = rep(100, times = nrow(T.AGI)), 
                          annotation.1 = rep(100, times = nrow(T.AGI)), annotation.2 = rep(100, times = nrow(T.AGI)))
i <- 1
for(i in i:nrow(T.AGI)){
  T.AGI$annotation.1[i] <- TAIR10$annotation[T.AGI$AGI.1[i] == TAIR10$AGI]
  T.AGI$annotation.2[i] <- TAIR10$annotation[T.AGI$AGI.2[i] == TAIR10$AGI]
  i <- i+1
}
i <- 1
for(i in i:nrow(T.data)){
  T.AGI$PeakIndex[T.AGI$AGI.1 == T.data$sAGI[i] & T.AGI$AGI.2 == T.data$tAGI[i]] <- T.data$PeakIndex[i]
  T.AGI$PeakIndex[T.AGI$AGI.2 == T.data$sAGI[i] & T.AGI$AGI.1 == T.data$tAGI[i]] <- T.data$PeakIndex[i]
  T.AGI$PeakMax[T.AGI$AGI.1 == T.data$sAGI[i] & T.AGI$AGI.2 == T.data$tAGI[i]] <- T.data$PeakMax[i]
  T.AGI$PeakMax[T.AGI$AGI.2 == T.data$sAGI[i] & T.AGI$AGI.1 == T.data$tAGI[i]] <- T.data$PeakMax[i]
  T.AGI$sameTFfamily[T.AGI$AGI.1 == T.data$sAGI[i] & T.AGI$AGI.2 == T.data$tAGI[i]] <- T.data$diffTFfamily[i]
  T.AGI$sameTFfamily[T.AGI$AGI.2 == T.data$sAGI[i] & T.AGI$AGI.1 == T.data$tAGI[i]] <- T.data$diffTFfamily[i]
}
T.matrix <- T.AGI[, c(3, 6, 7)] %>% spread(key = annotation.1, value = PeakIndex)
rownames(T.matrix) <- T.matrix$annotation.2
T.matrix <- T.matrix[, 2:ncol(T.matrix)]
res <- hclust(dist(T.matrix), method = "ward.D2")
T.AGI$annotation.1 <- factor(x = T.AGI$annotation.1,
                      levels = res$labels[res$order],
                      ordered = TRUE)
T.AGI$annotation.2 <- factor(x = T.AGI$annotation.2,
                      levels = res$labels[res$order],
                      ordered = TRUE)
T.AGI$PeakIndex[T.AGI$PeakIndex > 50] <- 50
heatmap.plot <- ggplot(NULL)
heatmap.plot <- heatmap.plot + theme_bw()
heatmap.plot <- heatmap.plot + geom_tile(data = T.AGI, aes(x = annotation.1, y = annotation.2, fill = PeakIndex), color = "white", size = 0.05)
heatmap.plot <- heatmap.plot + geom_tile(data = T.AGI %>% filter(PeakIndex > 2 & sameTFfamily == FALSE, PeakMax > 50), aes(x = annotation.1, y = annotation.2, fill = PeakIndex), color = "black", size = 0.3)
heatmap.plot <- heatmap.plot + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
heatmap.plot <- heatmap.plot + theme(axis.text = element_text(size = 8))
heatmap.plot + theme(legend.position = 'none')
#heatmap.plot <- heatmap.plot + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap.plot <- heatmap.plot + scale_fill_gradient(low = "white", high = "red", name="PeakIndex")
heatmap.plot <- heatmap.plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
heatmap.plot <- heatmap.plot + theme(axis.text = element_text(size = 10))
plot(heatmap.plot)
ggsave(filename = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/Image/UMAPCluster13_Fig/Heatmap.png", plot = heatmap.plot, width = 10, height = 6)