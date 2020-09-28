#"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/ATTEDII_EnrichemntOfTFsTargetGenes.R"
before <- proc.time()
install.packages("tidyr")
#package----
library(stringr)
library(dplyr)
library(tidyr)
#input data----
MasterTable <- readRDS("~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/20200311ATTEDII_MCLNumNodeTable.rds")
AtRegNet <- read.table("~/bigdata/yasue/AGRIS/AtRegNet/20190729AtRegNet_modified.txt", sep = "\t", quote = "", fill = T, stringsAsFactors = F, header = T)
AtRegNet <- AtRegNet[nchar(AtRegNet$TFLocus) == 9 & nchar(AtRegNet$TargetLocus) == 9, ]
TAIR10 <- readRDS("~/bigdata/yasue/TAIR10_ShortName.rds")
AGRIS <- read.table("~/bigdata/yasue/AGRIS/AGRIS_TFLIST/AGRIS_TFList_modified.txt", sep = "\t", quote = "", fill = T, stringsAsFactors = F)
AGRIS$V2 <- toupper(AGRIS$V2)
#processing data----
AtRegNet$TFLocus <- toupper(AtRegNet$TFLocus)
AtRegNet$TargetLocus <- toupper(AtRegNet$TargetLocus)
#Enrichment of TFs Target genes----
T.MCLNum <- max(as.numeric(MasterTable$MCLNum), na.rm = T)
TF.AGI <- unique(AtRegNet$TFLocus)
names(TF.AGI) <- TAIR10$annotation[match(TF.AGI, TAIR10$AGI)]
allgenes <- MasterTable$AGI
N <- length(allgenes)
T.data <- c()
j <- 1
for(j in j:length(TF.AGI)){
  T.target <- unique(AtRegNet$TargetLocus[AtRegNet$TFLocus == TF.AGI[j]])
  M <- length(T.target)
  T.pvalue <- c()
  k <- 1
  for(k in k:T.MCLNum){
    T.Node <- MasterTable %>% filter(MCLNum == k) %>% select("AGI") %>% unlist(use.names = F)
    n <- length(T.Node)
    x <- length(intersect(T.Node, T.target))
    T.pvalue <- c(T.pvalue, phyper(x-1, M, N-M, n, lower.tail = F))
    k <- k+1
  }
  T.data <- rbind(T.data, data.frame(MCLNum = 1:T.MCLNum,
                                     AGI = rep(TF.AGI[j], times = T.MCLNum),
                                     TF = rep(names(TF.AGI)[j], times = T.MCLNum),
                                     family = rep(AGRIS$V1[match(TF.AGI[j], AGRIS$V2)], times = T.MCLNum),
                                     pvalue = T.pvalue,
                                     qvalue = p.adjust(T.pvalue, method = "BH"),
                                     stringsAsFactors = F
  ))
  print(length(TF.AGI)-j)
  j <- j+1
}
#save data----
saveRDS(T.data, "~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_EnrichemntOfTFsTargetList.rds")
#elapsed time------------------------------------------
after <- proc.time()
print(after - before)#3240.220 sec, 54 min
#remove object----
rm(list = ls())