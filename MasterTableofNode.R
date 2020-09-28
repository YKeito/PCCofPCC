#"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/MasterTableofNode.R"
before <- proc.time()
#package----------------------------------------------------
library(dplyr)
library(stringr)
#input data-------------------------------------------------
CY.FPKM <- readRDS("~/Nakano_RNAseq/network_analysis/.RData/RDS/allCYFPKM.rds")
allgenes <- rownames(CY.FPKM)
NodeTable <- read.table("~/bigdata/yasue/PCCOfPCC/CY/Table/NetworkAnalyzerResults/CY15NodeTable_default.csv", sep = ",", header = T, stringsAsFactors = F)
CY.DEGs <- read.table("~/Nakano_RNAseq/network_analysis/base/CY_base/allRNASeq_union.txt", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
#Gene List of published Data
Botrytis_cinerea <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/defense/pathogen infection/old/B.cinerea.txt", sep = "\t", header = T, stringsAsFactors = F)
colnames(Botrytis_cinerea) <- "AGI"
PstDC3000_WT <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/defense/pathogen infection/pseudomonas syringae.txt", sep = "\t", header = T, stringsAsFactors = F)
colnames(PstDC3000_WT) <- "AGI"
PstDC3000_hrp <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/defense/pathogen infection/PstDC3000_hrp-_DEGs.txt", sep = "\t", header = T, stringsAsFactors = F)
MeJA_DEGs <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/MeJA_DEGs.txt", sep = "\t", header = T, stringsAsFactors = F)
SA_DEGs <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/SA_DEGs.txt", sep = "\t", header = T, stringsAsFactors = F)
SA_DEGs_1e3 <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/SA_DEGs_1e3.txt", sep = "\t", header = T, stringsAsFactors = F)
BTH <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/defense/BTH_affected_AGI.txt", sep = "\t", header = T, stringsAsFactors = F)
#processing data--------------------------------------------
#CY DEGs
T.CY <- c("CY15", "CY16", "CY20")
T.time <- c("1h", "3h", "12h", "24h")
allCY <- c()
obnames <- c()
i <- 1
for(i in i:length(T.CY)){
  n <- 1
  for(n in n:length(T.time)){
    temp <- CY.DEGs %>% select(starts_with(paste0(T.CY[i], "_",T.time[n]))) %>% select(ends_with("_q_value"))
    T.data <- rep("No", times = length(allgenes))
    names(T.data) <- allgenes
    T.data[match(rownames(temp)[temp < 0.05], names(T.data))] <- "Yes"
    allCY <- cbind(allCY, T.data)
    obnames <- c(obnames, paste0(T.CY[i], "_", T.time[n]))
    n <- n+1
  }
  i <- i+1
}
colnames(allCY) <- paste0(obnames, "_FDR0.05")
#Gene List
temp <- list(Botrytis_cinerea, PstDC3000_WT, PstDC3000_hrp, MeJA_DEGs, SA_DEGs, SA_DEGs_1e3, BTH)
obnames <- c("Botrytis_cinerea", "PstDC3000_WT", "PstDC3000_hrp", "MeJA_DEGs", "SA_DEGs", "SA_DEGs_1e3", "BTH")
allsample <- c()
i <- 1
for(i in i:length(obnames)){
  T.data <- rep("No", times = length(allgenes))
  names(T.data) <- allgenes
  T.data[match(temp[[i]][, "AGI"], names(T.data))] <- "Yes"
  allsample <- cbind(allsample, T.data)
  print(i)
  i <- i+1
}
colnames(allsample) <- obnames
#TF
TF_family <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/TF/arabidopsis_TF_family.txt", sep = "\t", header = T, stringsAsFactors = F)
T.data <- rep("No", times = length(allgenes))
names(T.data) <- allgenes
T.data[match(TF_family$AGI, names(T.data))] <- TF_family$TF
#MCL
MCLNum <- rep(NA, times = length(allgenes))
names(MCLNum) <- allgenes
MCLNum[match(NodeTable$name, names(MCLNum))] <- NodeTable$X__mclCluster
#Degree
Degree <- rep(0, times = length(allgenes))
names(Degree) <- allgenes
Degree[match(NodeTable$name, names(Degree))] <- NodeTable$Degree
#BetweennessCentrality
BC <- rep(0, times = length(allgenes))
names(BC) <- allgenes
BC[match(NodeTable$name, names(BC))] <- NodeTable$BetweennessCentrality
#MasterTable
MasterTable <- data.frame(AGI = allgenes,
                          TF = T.data,
                          MCLNum = MCLNum,
                          Degree = Degree,
                          BetweennessCentrality = BC,
                          allsample,
                          allCY,
                          stringsAsFactors = F
)
rownames(MasterTable) <- c()
#save data--------------------------------------------
saveRDS(MasterTable, "~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_MCLNumNodeTable.rds")
write.table(MasterTable, file = "~/bigdata/yasue/PCCOfPCC/CY/Table/MasterTableOfNode.txt", append=F, quote = F, sep = "\t", row.names = F)
#elapsed time----------------------------------------
after <- proc.time()
print(after - before)#1.427 sec
#remove object----------------------------------------
rm(list = ls())
