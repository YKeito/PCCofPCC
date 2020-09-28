"~/UMAP_script/05ATTEDII_UMAP_GO&TFtarget.R"
#package----
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(umap)
library(Rtsne)
library(tidyr)
#input data----
MasterTable <- readRDS(file = "~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_MCLNumNodeTable.rds")
ATTEDII_EnrichemntOfTFsTargetList <- readRDS("~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_EnrichemntOfTFsTargetList.rds")
T.ATTEDII_EnrichemntOfTFsTargetList <- ATTEDII_EnrichemntOfTFsTargetList %>% select(-pvalue, -AGI, -family)
colnames(T.ATTEDII_EnrichemntOfTFsTargetList) <- c("MCLNum", "Sample", "value")
allGOTable <- readRDS(file = "~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_allGOTable.rds")
allGOTable <-  allGOTable %>% select(-pvalue)
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
T.data <- rbind(allGOTable,
                T.ATTEDII_EnrichemntOfTFsTargetList
)
T.data$value <- -log10(T.data$value)
T.data <- T.data[T.data$MCLNum <= T.MCLNum, ]
T.data <- T.data %>% spread(key = Sample, value = value, convert = T)
rownames(T.data) <- formatC(T.data$MCLNum, width = "3", flag = "0")
T.data <- T.data %>% select(-MCLNum)
#いずれかのクラスターでターゲット遺伝子がエンリッチしているTFのリストアップ----
temp <- colnames(T.data)[apply(T.data >= -log10(0.05), MARGIN = 2, FUN = sum) != 0]
temp <- intersect(temp, ATTEDII_EnrichemntOfTFsTargetList$TF)
T.family <- ATTEDII_EnrichemntOfTFsTargetList$family[match(temp, ATTEDII_EnrichemntOfTFsTargetList$TF)]
names(T.family) <- ATTEDII_EnrichemntOfTFsTargetList$TF[match(temp, ATTEDII_EnrichemntOfTFsTargetList$TF)]
#T.data, T.family以外のobject削除----
rm(list = setdiff(ls(), c("T.data", "T.family")))
#umap----
#MCLNum----
T.config <- umap.defaults
T.config$random_state <- 0
T.config$n_neighbors <- 5
T.umap <- umap(T.data, config = T.config)
df.umap <- data.frame(x = T.umap$layout[, 1], 
                      y = T.umap$layout[, 2], 
                      Sample = rownames(T.data)
)
umap <- ggplot(df.umap, aes(x = x, y = y, color = Sample, label = Sample))
umap <- umap + geom_point()
umap <- umap + theme(legend.position = "top")
umap <- umap + theme(legend.title = element_blank())
umap <- umap + theme(legend.text = element_text(size=5))
umap <- umap + theme(legend.position = 'none')
umap <- umap + geom_text_repel(size = 2)
plot(umap)
#umap <- umap + geom_text(aes(label = rownames(T.data)), size = 3, hjust = 1)
ggsave("~/bigdata/yasue/PCCOfPCC/ATTEDII/Image/EnrichmentOfGenesList/Umap_MCLNum.png", umap)
#Sample----
T.data <- t(T.data)
T.config <- umap.defaults
T.config$random_state <- 0
T.config$n_neighbors <- 50
T.umap <- umap(T.data, config = T.config)
df.umap <- data.frame(x = T.umap$layout[, 1], 
                      y = T.umap$layout[, 2], 
                      Sample = rownames(T.data)
)
res <- dbscan(x = df.umap[, 1:2], eps = 0.28, minPts = 3)
df.umap <- df.umap %>% mutate(cluster = res$cluster)    # dataにクラスターの結果を追加
df.umap$cluster <- factor(df.umap$cluster)    # cluster列をfactor型に変換
umap <- ggplot(df.umap, aes(x = x, y = y, color = cluster, label = Sample))
umap <- umap + geom_point()
#umap <- umap + geom_text_repel(size = 1)
umap <- umap + geom_text(aes(label = rownames(T.data)), size = 2, vjust = 2.5)
#save data----
saveRDS(df.umap, "~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/UMAPClusterList.rds")
ggsave(paste0("~/bigdata/yasue/PCCOfPCC/ATTEDII/Image/EnrichmentOfGenesList/UMAP/Umap_Sample_n_neighbors", T.config$n_neighbors, "dbscan", res$eps, ".png"), umap, width = 8, height = 6)
write.table(x = df.umap, file = paste0("~/bigdata/yasue/PCCOfPCC/ATTEDII/Table/ATTEDII_", "Umap_Sample_n_neighbors", T.config$n_neighbors, "dbscan", res$eps, ".txt"), sep = "\t", quote = F, row.names = F)
#remove object----
rm(list = ls())