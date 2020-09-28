#"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/SummaryOfPoPResults"
before <- proc.time()
#package----
library(dplyr)
library(stringr)
#input data----
allPCCofPCCTable <- readRDS("~/bigdata/yasue/PCCOfPCC/CY/Table/PCCofPCC/naomit/allPCCofPCCTable.rds")
MasterTable <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_MCLNumNodeTable.rds")
allGOTable <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_allGOTable.rds")
MotifID_consensus <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CT15EnrichedMotif.rds")
CY.FPKM <- readRDS("~/Nakano_RNAseq/network_analysis/.RData/RDS/allCYFPKM.rds")
CY15_EnrichemntOfGenesList <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_EnrichemntOfGenesList.rds")
CY15_EnrichemntOfGenesList <- CY15_EnrichemntOfGenesList[CY15_EnrichemntOfGenesList$MCLNum < 115, ]
#processing data----
T.MCLNum <- rownames(allPCCofPCCTable)
T.MCLNum <- T.MCLNum[as.numeric(T.MCLNum) < 115]
CY15_EnrichemntOfGenesList <- CY15_EnrichemntOfGenesList %>% select("MCLNum", contains("qvalue:CY15"))
CY15_EnrichemntOfGenesList <- CY15_EnrichemntOfGenesList[match(T.MCLNum, CY15_EnrichemntOfGenesList$MCLNum), ]
CY15_EnrichemntOfGenesList <- CY15_EnrichemntOfGenesList[, 2:ncol(CY15_EnrichemntOfGenesList)] < 0.05
T.time <- c("1h", "3h", "12h", "24h")
T.CY <- paste0(rep(c("DMSO_", "CY15_"), each =4), T.time)
T.Node <- c()
T.condition <- c()
T.AvgExp <- c()
T.GO <- c()
T.AME <- c()
j <- 1
for(j in j:length(T.MCLNum)){
  T.data <- allGOTable %>% filter(MCLNum == T.MCLNum[j])
  T.data <- T.data$GOTerm[order(T.data$qvalue)][1:3]
  T.GO <- c(T.GO, paste(T.data, collapse = " | "))
  T.AGI <- MasterTable %>% filter(MCLNum == T.MCLNum[j]) %>% select("AGI") %>% unlist() %>% as.vector()
  T.average <- list()
  i <- 1
  for(i in i:length(T.CY)){
    T.average <- c(T.average, list(CY.FPKM[match(T.AGI, rownames(CY.FPKM)), ] %>% select(starts_with(T.CY[i])) %>% apply(MARGIN = 1, FUN = mean)))
    names(T.average)[[i]] <- T.CY[i]
    i <- i+1
  }
  T.average <- data.frame(T.average, stringsAsFactors = F)
  T.FPKM <- list()
  i <- 1
  for(i in i:length(T.time)){
    T.data <- T.average %>% select(ends_with(T.time[i]))
    T.FPKM <- c(T.FPKM, list(log2(c(T.data[, 2]+0.01)/c(T.data[, 1]+0.01))))
    names(T.FPKM[[i]]) <- rownames(T.average)
    i <- i+1
  }
  names(T.FPKM) <- paste0("CY15", T.time, "/DMSO", T.time)
  T.FPKM <- data.frame(T.FPKM, stringsAsFactors = F)
  T.FPKM <- T.FPKM %>% apply(MARGIN = 2, FUN = mean) %>% formatC(digits = 3)
  T.AvgExp <- rbind(T.AvgExp, T.FPKM)
  names(T.AvgExp)[j] <- T.MCLNum[j]
  T.Node <- c(T.Node, sum(MasterTable$MCLNum == T.MCLNum[j], na.rm = T))
  names(T.Node)[j] <- T.MCLNum[j]
  T.data <- allPCCofPCCTable[rownames(allPCCofPCCTable) == T.MCLNum[j], ] %>% select(ends_with("qvalue"))
  T.data <- colnames(T.data)[which(T.data < 0.05)]
  if(length(T.data) == 0){
    T.condition <- c(T.condition, NA)
    names(T.condition)[j] <- T.MCLNum[j]
  }else{
    T.condition <- c(T.condition, paste(str_split(T.data, pattern = "_qva", simplify = T)[, 1], collapse = "|"))
    names(T.condition)[j] <- T.MCLNum[j]
  }
  T.data <- MotifID_consensus$motif_alt_ID[MotifID_consensus$MCLNum == T.MCLNum[j]]
  if(length(T.data) !=0){
    T.AME <- c(T.AME, paste(T.data, collapse = " | "))
  }else{
    T.AME <- c(T.AME, NA)
  }
  print(length(T.MCLNum)-j)
  j <- j+1
}
colnames(T.AvgExp) <- paste0("Avglog2FC:CY15_", T.time)
colnames(CY15_EnrichemntOfGenesList) <- paste0("DEGs:CY15_", T.time)
T.data <- data.frame(MCLNum = T.MCLNum,
                     NumNodes = T.Node,
                     GOTerm = T.GO,
                     PoP_FDR005 = T.condition,
                     AME = T.AME,
                     CY15_DEGs = CY15_EnrichemntOfGenesList,
                     T.AvgExp,
                     stringsAsFactors = F
                     )
#save data----
saveRDS(T.data, file = "~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_Summary.rds")
write.table(T.data, file = "~/bigdata/yasue/PCCOfPCC/CY/Table/CY15_Summary.txt", sep = "\t", quote = F, row.names = F)
#elapsed time----
after <- proc.time()
print(after - before)#2.282 sec
#remove object----
rm(list = ls())
