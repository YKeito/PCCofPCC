#"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/ATTEDII_GenesListHeatmap.R"
#package----
library(ggplot2)
library(stringr)
library(gplots)
library(RColorBrewer)
#input data----
MasterTable <- readRDS(file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_MCLNumNodeTable.rds")
ATTEDII_EnrichemntOfTFsTargetList <- readRDS("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/EnrichmentAnalysis/ATTEDII_EnrichemntOfTFsTargetList.rds")
T.ATTEDII_EnrichemntOfTFsTargetList <- ATTEDII_EnrichemntOfTFsTargetList %>% select(-qvalue, -AGI, -family)
colnames(T.ATTEDII_EnrichemntOfTFsTargetList) <- c("MCLNum", "Sample", "value")
allGOTable <- readRDS(file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/EnrichmentAnalysis/ATTEDII_allGOTable.rds")
allGOTable <-  allGOTable %>% select(-qvalue)
colnames(allGOTable) <- c("MCLNum", "Sample", "value")
Hcluster.Results <- readRDS("~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_GenesListDendrogram.rds")
SHcluster.Results <- readRDS("~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_SampleDendrogram.rds")
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
T.data <- rbind(allGOTable %>% mutate(group = rep("GO", times = nrow(allGOTable))),
                T.ATTEDII_EnrichemntOfTFsTargetList %>% mutate(group = rep("TF Target", times = nrow(T.ATTEDII_EnrichemntOfTFsTargetList)))
)
T.data$value <- -log10(T.data$value)
T.data <- T.data[T.data$MCLNum <= T.MCLNum, ]
T.data$MCLNum <- formatC(T.data$MCLNum, width = "3", flag = "0")
#Heatmap----
#MCLNum
T.data$MCLNum <- factor(x = T.data$MCLNum,
                        levels = Hcluster.Results$MCLNum[Hcluster.Results$order],
                        ordered = TRUE)
T.data$Sample <- factor(x = T.data$Sample,
                        levels = as.factor(SHcluster.Results$Sample[SHcluster.Results$order]),
                        ordered = TRUE)
#Create heatmap plot
T.data$value[T.data$value >= 10] <- 10
heatmap.plot <- ggplot(data = T.data, aes(x = MCLNum, y = Sample, fill = value))
heatmap.plot <- heatmap.plot +geom_tile(color = "black", size = 0.01)
heatmap.plot <- heatmap.plot + scale_fill_gradient2(high="seagreen", low="white")
heatmap.plot <- heatmap.plot + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
heatmap.plot <- heatmap.plot + theme(axis.ticks.y = element_blank())
heatmap.plot <- heatmap.plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
heatmap.plot <- heatmap.plot + theme(legend.position = 'none')
heatmap.plot <- heatmap.plot + theme(axis.text=element_text(size=6))
heatmap.plot <- heatmap.plot + theme(axis.title.y = element_blank()) #軸名消す
heatmap.plot <- heatmap.plot + facet_grid(group~., scales="free_y", space="free", switch = "y")
heatmap.plot <- heatmap.plot + theme(strip.text.y = element_text(angle = 270))
heatmap.plot <- heatmap.plot + theme(strip.text = element_text(size = 12))
plot(heatmap.plot)
#save data----
ggsave("~/bigdata/yasue/PCCOfPCC/ATTEDII/Image/EnrichmentOfGenesList/ATTEDII_GenesListHeatmap.png", heatmap.plot, width = 12, height = 6)
#check:VennDiagram----
T.data <- ATTEDII_EnrichemntOfGenesList %>% spread(key = Sample, value = value, convert = T)
T.data <- T.data %>% filter(MCLNum <= T.MCLNum)
T.data <- T.data %>% select(-contains("BTH"), -contains("MCLNum"))
publish <- rownames(T.data)[which(apply(T.data %>% select(-contains("CY")) < 0.05, MARGIN = 1, FUN = sum) != 0)]
CY <- rownames(T.data)[which(apply(T.data %>% select(contains("CY")) < 0.05, MARGIN = 1, FUN = sum) != 0)]
CY15 <- rownames(T.data)[which(apply(T.data %>% select(contains("CY15")) < 0.05, MARGIN = 1, FUN = sum) != 0)]
CY16 <- rownames(T.data)[which(apply(T.data %>% select(contains("CY16")) < 0.05, MARGIN = 1, FUN = sum) != 0)]
CY20 <- rownames(T.data)[which(apply(T.data %>% select(contains("CY20")) < 0.05, MARGIN = 1, FUN = sum) != 0)]
T.CY <- list(CY15, CY16, CY20)
names(T.CY) <- c("CY15", "CY16", "CY20")
j <- 1
for(j in j:length(T.CY)){
  T.Venn <- list(publish, T.CY[[j]])
  names(T.Venn) <- c("published", names(T.CY)[j])
  ATTED.Venn <- venn(T.Venn)
  names(attr(ATTED.Venn,"intersections")) <- paste0(names(attr(ATTED.Venn,"intersections")), ",")
  T.data <- data.frame(CY = str_split(names(unlist(attr(ATTED.Venn,"intersections"))), pattern = ",", simplify = T)[, 1],
                       MCLNum = unlist(attr(ATTED.Venn,"intersections"), use.names = F),
                       stringsAsFactors = F)
  write.table(T.data, file = paste0("~/bigdata/yasue/PCCOfPCC/ATTEDII/Table/PublishedAnd", names(T.CY)[j], "VennDiagram.txt"), sep = "\t", quote = F, row.names = F)
  ggsave(filename = paste0("~/bigdata/yasue/PCCOfPCC/ATTEDII/Image/EnrichmentOfGenesList/PublishedAnd", names(T.CY)[j], "VennDiagram.png"), plot = plot(ATTED.Venn), width = 6, height = 4)
  j <- j+1
}
#remove object----
rm(list = ls())