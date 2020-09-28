#"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/ATTEEDII_MCLNumDendrogram.R"
#package----
library(reshape2)
library(factoextra)
library(ggplot2)
library(igraph)
library(dplyr)
library(stringr)
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
#must----
T.data <- rbind(allGOTable,
                T.ATTEDII_EnrichemntOfTFsTargetList
)
T.data$value <- -log10(T.data$value)
T.data <- T.data[T.data$MCLNum <= T.MCLNum, ]
T.data <- T.data %>% spread(key = Sample, value = value, convert = T)
rownames(T.data) <- formatC(T.data$MCLNum, width = "3", flag = "0")
T.data <- T.data %>% select(-MCLNum)
#Hclust----
#MCLNum
res.hc <- eclust(x = T.data, k = 10,
                 "hclust",
                 method = "euclidean",
                 graph = FALSE
)
#Sample
Sample.hc <- eclust(x = t(T.data), k = 20, 
                    "hclust",
                    method = "euclidean",
                    graph = FALSE
)
#Dendrogram----
#MCLNum
g <- fviz_dend(res.hc,
               cex = 0.5,
               color_labels_by_k = TRUE,
               show_labels = TRUE,
               ggtheme = theme_classic(),
               horiz = FALSE,
               #rect = TRUE,
               rect_fill = TRUE,
               type = "rectangle",
               main = NULL
)
g <- g + theme(axis.text=element_text(size=20))
g <- g + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
T.HCL <- data.frame(MCLNum = names(res.hc$cluster),
                    Hcluster = res.hc$cluster,
                    order = res.hc$order,
                    row.names = NULL)
#Sample
g.sample <- fviz_dend(Sample.hc,
                      cex = 0.5,
                      color_labels_by_k = TRUE,
                      show_labels = TRUE,
                      ggtheme = theme_classic(),
                      horiz = FALSE,
                      #rect = TRUE,
                      rect_fill = TRUE,
                      type = "rectangle",
                      main = NULL
)
g.sample <- g.sample + theme(axis.text=element_text(size=10))
g.sample <- g.sample + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
Sample.HCL <- data.frame(Sample = names(Sample.hc$cluster),
                         Hcluster = Sample.hc$cluster,
                         order = Sample.hc$order,
                         row.names = NULL)
#option----
#k:クラスター数の指定
#k_colors, palette:色指定
#show_labels:横軸の名前を表示するか
#color_labels_by_k:クラスターのグループに応じて色を付けるか
#horiz縦軸と横軸を反転するかしないか
#rect:クラスターを囲む
#rect_fill:rectで囲んだ範囲を色塗する
#type:図の形式の指定
#lwd:dendrogramの線の太さ変更
#label_cols:下の文字の色指定
#save----
#MCLNum
ggsave("~/bigdata/yasue/PCCOfPCC/ATTEDII/Image/EnrichmentOfGenesList/ATTEDII_MCLNumDendrogram.png", g, width = 20, height = 10)
saveRDS(T.HCL, "~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_GenesListDendrogram.rds")
write.table(T.HCL, "~/bigdata/yasue/PCCOfPCC/ATTEDII/Table/ATTEDII_GenesListDendrogram.txt", sep = "\t", append = F, quote = F, row.names = F, col.names = T)
#Sample
ggsave("~/bigdata/yasue/PCCOfPCC/ATTEDII/Image/EnrichmentOfGenesList/ATTEDII_SampleDendrogram.png", g.sample, width = 20, height = 10)
saveRDS(Sample.HCL, "~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_SampleDendrogram.rds")
write.table(Sample.HCL, "~/bigdata/yasue/PCCOfPCC/ATTEDII/Table/ATTEDII_SampleDendrogram.txt", sep = "\t", append = F, quote = F, row.names = F, col.names = T)