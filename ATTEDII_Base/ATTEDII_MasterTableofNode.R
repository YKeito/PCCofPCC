"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/ATTEDII_MasterTableofNode.R"
before <- proc.time()
#package----------------------------------------------------
library(dplyr)
library(stringr)
#input data-------------------------------------------------
NodeTable <- read.table("~/bigdata/yasue/ATTEDII/Microarray/NodeTable.csv", sep = ",", header = T, stringsAsFactors = F)
ATTEDII.AGI <- read.table("~/bigdata/yasue/ATTEDII/Microarray/Convert_EnterGeneID_to_TAIRID.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")
allgenes <- ATTEDII.AGI$To
allgenes <- sort(allgenes)
#Gene List of published Data
TF_family <- read.table("~/bigdata/yasue/AGRIS/AGRIS_TFLIST/AGRIS_TFList_modified.txt", sep = "\t", stringsAsFactors = F)
TF_family$V2 <- TF_family$V2 %>% toupper()
T.AGI <- intersect(TF_family$V2, allgenes)
TF_family <- TF_family[match(T.AGI, TF_family$V2), ]
#processing data--------------------------------------------
obnames <- c("X__mclCluster", "Degree", "BetweennessCentrality")
NetworkStructure <- c()
i <- 1
for(i in i:length(obnames)){
  T.data <- rep("No", times = length(allgenes))
  names(T.data) <- allgenes
  T.data[match(NodeTable$name, names(T.data))] <- NodeTable[, obnames[i]]
  NetworkStructure <- cbind(NetworkStructure, T.data)
  print(i)
  i <- i+1
}
colnames(NetworkStructure) <- c("MCLNum", "Degree", "BetweennessCentrality")
#TF
T.data <- rep("No", times = length(allgenes))
names(T.data) <- allgenes
T.data[match(TF_family$V2, names(T.data))] <- TF_family$V1

#MasterTable
MasterTable <- data.frame(AGI = allgenes,
                          TF = T.data,
                          NetworkStructure,
                          stringsAsFactors = F,
                          row.names = NULL
                          )
#save data--------------------------------------------
#saveRDS(MasterTable, "~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/20200311ATTEDII_MCLNumNodeTable.rds")
#write.table(MasterTable, file = "~/bigdata/yasue/PCCOfPCC/ATTEDII/Table/20200311ATTEDII_MasterTableOfNode.txt", append=F, quote = F, sep = "\t", row.names = F)
#elapsed time----------------------------------------
after <- proc.time()
print(after - before)#1.489 sec
#remove object----------------------------------------
rm(list = ls())