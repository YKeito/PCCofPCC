"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/ATTEDII_EnrichemntOfTFsTargetList.R"
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
ATTEDII_EnrichemntOfTFsTargetList <- ATTEDII_EnrichemntOfTFsTargetList %>% select(-pvalue, -AGI, -family)
colnames(ATTEDII_EnrichemntOfTFsTargetList) <- c("MCLNum", "Sample", "value")
allGOTable <- readRDS(file = "~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_allGOTable.rds")
allGOTable <-  allGOTable %>% select(-pvalue)
colnames(allGOTable) <- c("MCLNum", "Sample", "value")
ATTEDII_EnrichemntOfGenesList <- readRDS("~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_EnrichemntOfGenesList.rds")
ATTEDII_EnrichemntOfGenesList <- ATTEDII_EnrichemntOfGenesList %>% select(-"pvalue")
ATTEDII_EnrichemntOfGenesList <- ATTEDII_EnrichemntOfGenesList[ATTEDII_EnrichemntOfGenesList$Sample != "BTH", ]
colnames(ATTEDII_EnrichemntOfGenesList) <- c("MCLNum", "Sample", "value")
#allAverageExp <- readRDS("~/bigdata/yasue/GEO/RDS/AverageExp/allAverageExp_MCLNum.rds")
#colnames(allAverageExp) <- c("MCLNum", "Sample", "value")
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
#must----
T.data <- rbind(ATTEDII_EnrichemntOfGenesList,
                allGOTable,
                ATTEDII_EnrichemntOfTFsTargetList
)
T.data$value <- -log10(T.data$value)
#T.data <- rbind(T.data, allAverageExp)
T.data <- T.data[T.data$MCLNum <= T.MCLNum, ]
T.data <- T.data %>% spread(key = Sample, value = value, convert = T)
rownames(T.data) <- formatC(T.data$MCLNum, width = "3", flag = "0")
T.data <- T.data %>% select(-MCLNum)
#umap----
#MCLNum----
T.config <- umap.defaults
T.config$random_state <- 0
T.config$n_neighbors <- 5
T.config$n_epochs <- 20
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
ggsave("~/bigdata/yasue/PCCOfPCC/ATTEDII/Image/EnrichmentOfGenesList/UMAP/Umap_MCLNum.png", umap)
#Sample----
T.data <- t(T.data)
T.config <- umap.defaults
T.config$random_state <- 0
#T.config$n_neighbors <- 10
#T.config$spread <- 20
#T.config$alpha <- 5
T.umap <- umap(T.data, config = T.config)
df.umap <- data.frame(x = T.umap$layout[, 1], 
                      y = T.umap$layout[, 2], 
                      Sample = rownames(T.data))
umap <- ggplot(df.umap, aes(x = x, y = y, color = Sample, label = Sample))
umap <- umap + geom_point()
#umap <- umap + geom_text_repel(size = 4)
umap <- umap + geom_text(aes(label = rownames(T.data)), size = 1, vjust = 2.5)
umap <- umap + theme(legend.title = element_blank())
umap <- umap + theme(legend.position = 'none')
#plot(umap)
ggsave("~/bigdata/yasue/PCCOfPCC/ATTEDII/Image/EnrichmentOfGenesList/UMAP/Umap_Sample.png", umap)