#"/home/yasue/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/FioreDB/UMAPCluster_FioreDB_Heatmap.R"
#load package----
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
#Input data----
#Individual project
I.data <- read.table("/home/yasue/bigdata/yasue/FioreDB/IndividualProject.txt", sep = "\t", header = T, stringsAsFactors = F)
I.data$Gene <- toupper(I.data$Gene)
I.data <- I.data %>% filter(Num.of.lines.showing.the.phenotype > 1)
#Bulk project
B.data <- read.table("/home/yasue/bigdata/yasue/FioreDB/BulkProject.txt", sep = "\t", header = T, stringsAsFactors = F)
B.data[is.na(B.data$Phenotype.code), 4:7] <- B.data[is.na(B.data$Phenotype.code), 5:8]
B.data <- B.data %>% select(-X)
B.data$Causal.gene <- toupper(B.data$Causal.gene)
B.data <- B.data[grep("AT", B.data$Causal.gene), ]
PhenotypeData <- data.frame(AGI = c(I.data$Gene, B.data$Causal.gene),
                            PhenotypeCode = c(I.data$Phenotype.code, B.data$Phenotype.code),
                            Tissue = c(I.data$Tissue, B.data$Tissue),
                            Phenotype = c(I.data$Phenotype, B.data$Phenotype),
                            stringsAsFactors = F
)
#UMAPClusterのTFリスト-----
DAPSeq.TFName <- readRDS("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/DAPSeq.TFName.rds")
T.AGI <- DAPSeq.TFName$AGI
T.ID <- unique(PhenotypeData$PhenotypeCode) %>% sort()
T.matrix <- data.frame(AGI = T.AGI,
                       UMAPCluster = DAPSeq.TFName$UMAPCluster,
                       temp = matrix(rep(0, times = length(T.AGI)*length(T.ID)), nrow = length(T.AGI), ncol = length(T.ID)),
                       stringsAsFactors = F)
colnames(T.matrix)[3:ncol(T.matrix)] <- paste0("PhenotypeID_", T.ID)
h <- 1
for(h in h:length(T.ID)){
  T.AGI <- PhenotypeData$AGI[PhenotypeData$PhenotypeCode == T.ID[h]] %>% unique()
  T.AGI <- intersect(T.AGI, T.data$AGI)
  if(length(T.AGI) != 0){
    T.matrix[match(T.AGI, T.matrix$AGI), 2+h] <- 1
  }
  print(h)
  h <- h+1
}
T.data <- T.matrix %>% gather(key = FioreID, value = value, PhenotypeID_0:PhenotypeID_999)
T.data <- T.data %>% mutate(Tissue = rep(0, times = nrow(T.data)),
                            phenotype = rep(0, times = nrow(T.data))
)
T.key <- paste0("PhenotypeID_", PhenotypeData$PhenotypeCode) %>% unique() %>% sort()
h <- 1
for(h in h:length(T.key)){
  T.PhenotypeData <- PhenotypeData %>% filter(PhenotypeCode == str_split(T.key[h], pattern = "_", simplify = T)[, 2]) %>% distinct()
  T.data$Tissue[grep(T.key[h], T.data$FioreID)] <- T.PhenotypeData$Tissue[1]
  T.data$phenotype[grep(T.key[h], T.data$FioreID)] <- T.PhenotypeData$Phenotype[1]
  h <- h+1
}
#HCL clustering----
res.matrix <- T.matrix[3:ncol(T.matrix)]
rownames(res.matrix) <- T.matrix$AGI
res.1 <- hclust(dist(res.matrix), method = "ward.D2")
T.data$AGI <- factor(x = T.data$AGI,
                     levels = res.1$labels[res.1$order],
                     ordered = TRUE)
res.2 <- hclust(dist(t(res.matrix)), method = "ward.D2")
T.data$FioreID <- factor(x = T.data$FioreID,
                         levels = res.2$labels[res.2$order],
                         ordered = TRUE)
T.data$UMAPCluster <- factor(T.data$UMAPCluster,
                             0:17)
#Heatmap----
g <- ggplot(T.data, aes(x = phenotype, y = AGI, fill = value))
g <- g + geom_tile(color = "black", size = 0.01)
g <- g + theme_bw()
g <- g + scale_fill_gradient2(limits = c(0, 1), low = "blue", high = "red", na.value = "white")
g <- g + facet_grid(UMAPCluster~Tissue, scales="free", space="free", switch = "y")
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g <- g + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
g <- g + theme(axis.ticks.y = element_blank())
g <- g + theme(axis.text.y = element_blank())
g <- g + theme(legend.position = 'none')
plot(g)
ggsave("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/Image/FioreDB_Phenotype_Heatmap.png", g, width = 12, height = 8)
#SortMatrix HCL clustering----
Uni.PhenotypeData <- PhenotypeData %>% arrange(PhenotypeCode) %>% distinct(PhenotypeCode, .keep_all = T)
T.Tissue <- c("", "", Uni.PhenotypeData$Tissue)
T.Phenotype <- c("", "", Uni.PhenotypeData$Phenotype)
T.matrix <- rbind(T.Tissue, T.Phenotype, T.matrix)
sort.matrix <- T.matrix[c(1:2, res.1$order+2), c(1:2, res.2$order+2)]
sort.matrix$UMAPCluster[3:nrow(sort.matrix)] <- paste0("UMAPCluster", sort.matrix$UMAPCluster[3:nrow(sort.matrix)])
saveRDS(object = sort.matrix, file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/FioreDB/UMAPCluster_FioreDB_SortMatrix.rds")
write.table(sort.matrix, file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/Table/FioreDB/UMAPCluster_FioreDB_SortMatrix.txt", sep = "\t", quote = F, row.names = F)