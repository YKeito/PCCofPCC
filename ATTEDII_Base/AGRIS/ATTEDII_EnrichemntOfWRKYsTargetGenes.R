"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/ATTEDII_EnrichemntOfWRKYsTargetGenes.R"
#package----
library(readr)
library(stringr)
library(dplyr)
#input data----
MasterTable <- readRDS("~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_MCLNumNodeTable.rds")
AgriNet <- read.table("~/bigdata/yasue/OpenData/AtRegNet_processed.csv", sep = ",", header = T, , quote = "", fill = T, stringsAsFactors = F)
#processing data----
#AgriNet <- read_delim("~/bigdata/yasue/OpenData/AtRegNet.csv", delim = ",", quote = "")
AgriNet$TFLocus <- toupper(AgriNet$TFLocus)
AgriNet$TFLocus <- gsub(AgriNet$TFLocus, pattern = " ", replacement = "")
AgriNet$TFName <- toupper(AgriNet$TFName)
AgriNet$TargetLocus <- toupper(AgriNet$TargetLocus)
AgriNet$TargetLocus <- gsub(AgriNet$TargetLocus, pattern = " ", replacement = "")
AgriNet$TargetName <- toupper(AgriNet$TargetName)
#select WRKY----
T.WRKY <- AgriNet[grep("WRKY", AgriNet$TFFamily), ]
T.WRKY$TFName <- gsub(T.WRKY$TFName, pattern=" ", replacement="")
T.WRKY$TFName <- gsub(T.WRKY$TFName, pattern="AT", replacement="")
T.AGI <- unique(T.WRKY$TFLocus)
WRKY.Target <- c()
i <- 1
for(i in i:length(T.AGI)){
  WRKY.Target <- rbind(WRKY.Target, T.WRKY %>% filter(TFLocus == T.AGI[i]) %>% select(TFName, TFLocus, TargetName, TargetLocus))
  i <- i+1
}
#Node ≧ 10----
T.MCLNum <- unique(MasterTable$MCLNum)
T.MCLNum <- T.MCLNum[!is.na(T.MCLNum)]
T.Node <- c()
j <- 1
for(j in j:length(T.MCLNum)){
  T.Node <- c(T.Node, MasterTable %>% filter(MCLNum == j) %>% select("AGI") %>% nrow())
  j <- j+1
}
T.MCLNum <- max(which(T.Node == 10))

#Enrichment of WRKYs Target genes----
allgenes <- MasterTable$AGI
N <- length(allgenes)
T.data <- c()
T.sample <- unique(WRKY.Target$TFName)
j <- 1
for(j in j:length(T.sample)){
  T.target <- WRKY.Target %>% filter(TFName == T.sample[j]) %>% select(TargetLocus) %>% unlist(use.names = F)
  M <- length(T.target)
  T.pvalue <- c()
  k <- 1
  for(k in k:T.MCLNum){
    T.Node <- MasterTable %>% filter(MCLNum == k) %>% select("AGI")
    n <- nrow(T.Node)
    x <- length(intersect(T.Node$AGI, T.target))
    T.pvalue <- c(T.pvalue, phyper(x-1, M, N-M, n, lower.tail = F))
    k <- k+1
  }
  T.data <- rbind(T.data, data.frame(MCLNum = 1:max(T.MCLNum),
                       WRKY = rep(T.sample[j], times = T.MCLNum),
                       T.pvalue,
                       p.adjust(T.pvalue, method = "BH"),
                       stringsAsFactors = F
                       ))
  print(length(T.sample)-j)
  j <- j+1
}
saveRDS(T.data, "~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_EnrichemntOfWRKYsTargetList.rds")
#remove object----
rm(list = ls())