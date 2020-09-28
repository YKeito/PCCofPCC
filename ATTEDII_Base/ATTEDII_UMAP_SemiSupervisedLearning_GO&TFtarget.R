"/home/yasue/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/ATTEDII_UMAP_SemiSupervisedLearning_GO&TFtarget.R"
#package----
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(umap)
library(dbscan)
library(tidyr)
#input data----
MasterTable <- readRDS(file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_MCLNumNodeTable.rds")
ATTEDII_EnrichemntOfTFsTargetList <- readRDS("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/EnrichmentAnalysis/ATTEDII_EnrichemntOfTFsTargetList.rds")
T.ATTEDII_EnrichemntOfTFsTargetList <- ATTEDII_EnrichemntOfTFsTargetList %>% select(-qvalue, -AGI, -family)
colnames(T.ATTEDII_EnrichemntOfTFsTargetList) <- c("MCLNum", "Sample", "value")
allGOTable <- readRDS(file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/EnrichmentAnalysis/ATTEDII_allGOTable.rds")
allGOTable <-  allGOTable %>% select(-qvalue)
colnames(allGOTable) <- c("MCLNum", "Sample", "value")
#processing data----
#Node ≧ 8----
T.MCLNum <- unique(MasterTable$MCLNum)
T.MCLNum <- T.MCLNum[!is.na(T.MCLNum)]
T.Node <- c()
j <- 1
for(j in j:length(T.MCLNum)){
  T.Node <- c(T.Node, MasterTable %>% filter(MCLNum == j) %>% select("AGI") %>% nrow())
  j <- j+1
}
T.MCLNum <- max(which(T.Node == 8))
#複数のデータをrbindで成型----
USL.data <- rbind(allGOTable,
                  T.ATTEDII_EnrichemntOfTFsTargetList
)
USL.data$value <- -log10(USL.data$value)
USL.data <- USL.data %>% filter(MCLNum <= T.MCLNum)
USL.data <- USL.data %>% spread(key = Sample, value = value, convert = T)
USL.data <- USL.data %>% select(-MCLNum)
#UMAP----
#MCLNum Unsupervised learning----
T.config <- umap.defaults
T.config$random_state <- 0
T.config$n_neighbors <- 8
umap.USL <- umap(USL.data, config = T.config)
df.umap <- data.frame(x = umap.USL$layout[, 1], 
                      y = umap.USL$layout[, 2], 
                      Sample = rownames(USL.data)
)
res <- dbscan(x = df.umap[, 1:2], eps = 0.3, minPts = 3)
df.umap <- df.umap %>% mutate(cluster = res$cluster)    # dataにクラスターの結果を追加
df.umap$cluster <- factor(df.umap$cluster)    # cluster列をfactor型に変換
umap <- ggplot(df.umap, aes(x = x, y = y, color = cluster, label = Sample))
umap <- umap + geom_point()
umap <- umap + ggtitle(paste0("n_neighbors", T.config$n_neighbors, "DBSCAN:eps = ", res$eps))
umap <- umap + theme(legend.position = 'none')
plot(umap)
#ggsave("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/Image/EnrichmentOfGenesList/20190809Umap_MCLNum.png", umap)
#Sample----
USL.data <- t(USL.data)
T.config <- umap.defaults
T.config$random_state <- 0
T.config$n_neighbors <- 40
umap.USL <- umap(USL.data, config = T.config)
df.umap <- data.frame(x = umap.USL$layout[, 1], 
                      y = umap.USL$layout[, 2], 
                      Sample = rownames(USL.data)
)
res <- dbscan(x = df.umap[, 1:2], eps = 0.4, minPts = 3)
df.umap <- df.umap %>% mutate(cluster = res$cluster)    # dataにクラスターの結果を追加
df.umap$cluster <- factor(df.umap$cluster)    # cluster列をfactor型に変換
umap <- ggplot(data = NULL)
umap <- umap + geom_point(data = df.umap, aes(x = x, y = y, color = cluster, label = Sample))
#umap <- umap + geom_point(data = df.umap, aes(x = x, y = y, label = Sample), color = "black")
umap <- umap + ggtitle(paste0("n_neighbors", T.config$n_neighbors, "DBSCAN:eps = ", res$eps))
umap <- umap + xlab("UMAP1") + ylab("UMAP2")
umap <- umap + theme_bw()
#umap <- umap + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) #軸名消す
#umap <- umap + theme(axis.text = element_blank()) #軸の値消す
#umap <- umap + theme(legend.position = 'none') #凡例消す
plot(umap)
#save data----
#ggsave(paste0("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/Image/EnrichmentOfGenesList/UMAP/FigUmap_Sample_n_neighbors", T.config$n_neighbors, "dbscan", res$eps, ".png"), umap, width = 8, height = 6)
#saveRDS(df.umap, "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/20190809UMAPClusterList.rds")
#write.table(x = df.umap, file = paste0("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/Table/20190809ATTEDII_", "Umap_Sample_n_neighbors", T.config$n_neighbors, "dbscan", res$eps, ".txt"), sep = "\t", quote = F, row.names = F)
#MCLNum Semi-Supervised learning----
#input data----
ATTEDII_EnrichemntOfGenesList <- readRDS("~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/EnrichmentAnalysis/ATTEDII_EnrichemntOfGenesList.rds")
ATTEDII_EnrichemntOfGenesList <- ATTEDII_EnrichemntOfGenesList %>% select(-"qvalue")
colnames(ATTEDII_EnrichemntOfGenesList) <- c("MCLNum", "Sample", "value")
ATTEDII_EnrichemntOfGenesList<- ATTEDII_EnrichemntOfGenesList %>% filter(MCLNum <= T.MCLNum)
ATTEDII_EnrichemntOfGenesList$value <- -log10(ATTEDII_EnrichemntOfGenesList$value)
ATTEDII_EnrichemntOfGenesList <- ATTEDII_EnrichemntOfGenesList %>% spread(key = Sample, value = value, convert = T)
#T.data <- ATTEDII_EnrichemntOfGenesList %>% select(-contains("CY"), -"MCLNum", -"SA_DEGs", -"PstDC3000_hrp") %>% t()
#rownames(T.data) <- paste0(0, 1:nrow(T.data), rownames(T.data))
ATTEDII_EnrichemntOfGenesList <- ATTEDII_EnrichemntOfGenesList %>% select(contains("CY")) #ATTEDII_EnrichemntOfGenesList <- ATTEDII_EnrichemntOfGenesList %>% select(contains("CY"), "BTH")
colnames(ATTEDII_EnrichemntOfGenesList) <- str_split(colnames(ATTEDII_EnrichemntOfGenesList), pattern = "_FDR", simplify = T)[, 1]
ATTEDII_EnrichemntOfGenesList <- ATTEDII_EnrichemntOfGenesList %>% t()
T.CY <- str_split(rownames(ATTEDII_EnrichemntOfGenesList), pattern = "_", simplify = T)
temp <- as.numeric(str_split(T.CY[, 2], pattern = "h", simplify = T)[, 1])
rownames(ATTEDII_EnrichemntOfGenesList) <- paste0(T.CY[, 1], "_", formatC(temp, width = 2, flag = "0"), "h")
ATTEDII_EnrichemntOfGenesList <- ATTEDII_EnrichemntOfGenesList[order(rownames(ATTEDII_EnrichemntOfGenesList)), ]
#ATTEDII_EnrichemntOfGenesList <- rbind(T.data, ATTEDII_EnrichemntOfGenesList)
#Data of Semisupervised learning----
df2.umap <- c()
h <- 1
for(h in h:nrow(ATTEDII_EnrichemntOfGenesList)){
  SSL.data <- matrix(data = ATTEDII_EnrichemntOfGenesList[h, ], ncol = ncol(ATTEDII_EnrichemntOfGenesList))
  umap.SSL <- predict(umap.USL, SSL.data)
  df2.umap <- rbind(df2.umap, data.frame(x = umap.SSL[, 1], 
                                         y = umap.SSL[, 2], 
                                         Sample = rownames(ATTEDII_EnrichemntOfGenesList)[h]
  ))
  h <- h+1
}
df2.umap$Sample <- factor(df2.umap$Sample, df2.umap$Sample[c(17, 15, 16, 1, 2, 4, 6, 3, 5, c(c(4, 6, 3, 5)+4), c(c(4, 6, 3, 5)+8))])
#factor(df2.umap$Sample, df2.umap$Sample[c(4, 2, 3, 1)]) 
#factor(df2.umap$Sample, df2.umap$Sample[c(2, 5, 3, 4, 1)]), 
#factor(df2.umap$Sample, df2.umap$Sample[c(2, 5, 3, 4, 1, 7, 9, 6, 8, 11, 13, 10, 12, 15, 17, 14, 16)])
umap2 <- umap + geom_point(data = df2.umap, aes(x = x, y = y, shape = Sample), size = 5)
umap2 <- umap2 + scale_shape_manual(values=1:nlevels(df2.umap$Sample))
umap2 <- umap2 + guides(shape = guide_legend(order = 2), colour = guide_legend(order = 1))
umap2 <- umap2 + theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8),
                       legend.key.size = unit(0.5, "cm"), legend.key.width = unit(0.5,"cm")
)
plot(umap2)
ggsave(paste0("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/Image/EnrichmentOfGenesList/UMAP/Fig2CY&PubSSLUmap_Sample_n_neighbors", T.config$n_neighbors, "dbscan", res$eps, ".png"), umap2, width = 8, height = 6)
#remove object----
rm(list = ls())