#"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/CY15MCLNum_GOanalysis.R"
before <- proc.time()
#package----
library(dplyr)
#input data------------------------------------------------------------------------------------------------------------------------------------
MasterTable <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_MCLNumNodeTable.rds")
GO_List <- readRDS("~/bigdata/yasue/GO_analysis/GO_List.rds")
ATH_GOSLim <- readRDS("~/bigdata/yasue/GO_analysis/ATH_GO_GOSLIM.rds")
#processing data-------------------------------------------------------------------------------------------------------------------------------
T.GOList <- unlist(GO_List)
N <- length(T.GOList)
UniGO <- unique(T.GOList)
T.MCLNum <- MasterTable$MCLNum
T.MCLNum <- T.MCLNum[!is.na(T.MCLNum)]
T.MCLNum <- T.MCLNum[as.numeric(T.MCLNum) < 115]
T.MCLNum <- unique(T.MCLNum)
allGOTable <- c()
j <- 1
for(j in j:length(T.MCLNum)){
  T.Node <- MasterTable %>% filter(MCLNum == T.MCLNum[j]) %>% select("AGI")
  T.subGO <- GO_List[match(T.Node$AGI, names(GO_List))]
  n <- length(unlist(T.subGO))
  T.pvalue <- c()
  T.AGI <- c()
  i <- 1
  for(i in i:length(UniGO)){
    x <- sum(unlist(T.subGO) == UniGO[i])
    T.AGI <- c(T.AGI, paste(names(T.subGO)[grep(UniGO[i], T.subGO)], collapse = "|"))
    M <- sum(UniGO[i] == T.GOList)
    T.pvalue <- c(T.pvalue, phyper(x-1, M, N-M, n, lower.tail = F))
    names(T.pvalue)[i] <- UniGO[i]
    i <- i+1
  }
  T.data <- data.frame(MCLNum = rep(T.MCLNum[j], times = length(T.AGI)),
                       GOID = UniGO,
                       GOTerm = ATH_GOSLim$`GO term`[match(UniGO, ATH_GOSLim$`GO ID`)],
                       pvalue = T.pvalue,
                       qvalue = p.adjust(T.pvalue, method = "BH"),
                       AGI = T.AGI,
                       stringsAsFactors = F
                       )
  T.data <- T.data[T.data$qvalue < 0.05, ]
  if(nrow(T.data) != 0){
    allGOTable <- rbind(allGOTable, T.data)
    title <- paste0("~/bigdata/yasue/PCCOfPCC/CY/Table/GOAnalysis/CY15_MCLNum", T.MCLNum[j], "_GOTable.txt")
    write.table(T.data, file = title, sep = "\t", quote = F, row.names = F)
  }
  print(length(T.MCLNum)-j)
  j <- j+1
}
#save data----
saveRDS(allGOTable, file = "~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_allGOTable.rds")
write.table(allGOTable, file = "~/bigdata/yasue/PCCOfPCC/CY/Table/GOAnalysis/CY15_allGOTable.txt", sep = "\t", quote = F, row.names = F)
#elapsed time------------------------------------------
after <- proc.time()
print(after - before)#2832.952 sec, 47.21587 min
#remove object-----------------------------------------
rm(list = ls())
