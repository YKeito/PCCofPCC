"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/ATTEDII_EnrichemntOfGenesList.R"
library(dplyr)
before <- proc.time()
#input data---------------------------------
MasterTable <- readRDS("~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/20200311ATTEDII_MCLNumNodeTable.rds")
#processing data----
T.MCLNum <- max(as.numeric(MasterTable$MCLNum), na.rm = T)
allgenes <- MasterTable$AGI
N <- length(allgenes)
T.data <- c()
T.sample <- colnames(MasterTable)
T.sample <- setdiff(T.sample, c("TF", "MCLNum", "Degree", "BetweennessCentrality"))
j <- 2
for(j in j:length(T.sample)){
  T.DEGs <- MasterTable %>% select("AGI", T.sample[j])
  T.DEGs <- T.DEGs$AGI[T.DEGs[, 2] == "Yes"]
  M <- length(T.DEGs)
  T.pvalue <- c()
  k <- 1
  for(k in k:T.MCLNum){
    T.Node <- MasterTable %>% filter(MCLNum == k) %>% select("AGI")
    n <- nrow(T.Node)
    x <- length(intersect(T.Node$AGI, T.DEGs))
    T.pvalue <- c(T.pvalue, phyper(x-1, M, N-M, n, lower.tail = F))
    k <- k+1
  }
  T.data <- rbind(T.data, data.frame(MCLNum = 1:T.MCLNum,
                                     Sample = rep(T.sample[j], times = T.MCLNum),
                                     pvalue = T.pvalue,
                                     qvalue = p.adjust(T.pvalue, method = "BH"),
                                     stringsAsFactors = F)
  )
  print(length(T.sample)-j)
  j <- j+1
}
#save data---------------------------------------------
saveRDS(T.data, "~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_EnrichemntOfGenesList.rds")
write.table(T.data, file = "~/bigdata/yasue/PCCOfPCC/ATTEDII/Table/ATTEDII_EnrichemntOfGenesList.txt", sep = "\t", quote = F, row.names = F)
#elapsed time------------------------------------------
after <- proc.time()
print(after - before)#105.467 sec
#remove object-----------------------------------------
rm(list = ls())