#"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/ATTEDIIMCLNum_GOAnalysis.R"
before <- proc.time()
#package----
library(dplyr)
#input data------------------------------------------------------------------------------------------------------------------------------------
MasterTable <- readRDS("~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_MCLNumNodeTable.rds")
GOSLIMTermProcess <- readRDS("~/bigdata/yasue/GO_analysis/GOList_GOSLIMTermProcess.rds")
ATHGOSLim.Process <- readRDS("~/bigdata/yasue/GO_analysis/ATH_GO_GOSLIM_Process.rds")
TAIR10.AGI <- readRDS("~/bigdata/yasue/TAIR10_AGI.rds")
#processing data-------------------------------------------------------------------------------------------------------------------------------
N <- length(TAIR10.AGI)
UniGO <- names(GOSLIMTermProcess)
T.MCLNum <- max(as.numeric(MasterTable$MCLNum), na.rm = T)
allGOTable <- c()
j <- 1
for(j in j:length(UniGO)){
  GO.AGI <- GOSLIMTermProcess[[which(names(GOSLIMTermProcess) == UniGO[j])]]
  allAGI.Process <- intersect(TAIR10.AGI, GO.AGI)
  M <- length(allAGI.Process)
  T.pvalue <- c()
  T.AGI <- c()
  i <- 1
  for(i in i:T.MCLNum){
    T.Node <- MasterTable %>% filter(MCLNum == i) %>% select(AGI) %>% unlist(use.names = F)
    n <- length(T.Node)
    SubAndGO.AGI <- intersect(T.Node, GO.AGI)
    x <- length(intersect(T.Node, SubAndGO.AGI))
    T.AGI <- c(T.AGI, paste(SubAndGO.AGI, collapse = "|"))
    T.pvalue <- c(T.pvalue, phyper(x-1, M, N-M, n, lower.tail = F))
    i <- i+1
  }
  T.data <- data.frame(MCLNum = 1:T.MCLNum,
                       Sample = rep(UniGO[j], times = T.MCLNum),
                       pvalue = T.pvalue,
                       qvalue = p.adjust(T.pvalue, method = "BH"),
                       stringsAsFactors = F
  )
  allGOTable <- rbind(allGOTable, T.data)
  print(length(UniGO)-j)
  j <- j+1
}
#save data----
saveRDS(allGOTable, file = "~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_allGOTable.rds")
write.table(allGOTable, file = "~/bigdata/yasue/PCCOfPCC/ATTEDII/Table/GOResults/ATTEDII_allGOTable.txt", sep = "\t", quote = F, row.names = F)
#elapsed time------------------------------------------
after <- proc.time()
print(after - before)#267.159 sec
#remove object-----------------------------------------
rm(list = ls())