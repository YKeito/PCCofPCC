#"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/EnrichemntOfGenesList.R"
library(dplyr)
before <- proc.time()
#input data---------------------------------
MasterTable <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_MCLNumNodeTable.rds")
#processing data----------------------------
T.MCLNum <- unique(MasterTable$MCLNum)
T.MCLNum <- T.MCLNum[!is.na(T.MCLNum)]
T.Node <- c()
j <- 1
for(j in j:length(T.MCLNum)){
  T.Node <- c(T.Node, MasterTable %>% filter(MCLNum == j) %>% select("AGI") %>% nrow())
  j <- j+1
}
T.MCLNum <- max(which(T.Node == 10))

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
  for(k in k:max(T.MCLNum)){
    T.Node <- MasterTable %>% filter(MCLNum == k) %>% select("AGI")
    n <- nrow(T.Node)
    x <- length(intersect(T.Node$AGI, T.DEGs))
    T.pvalue <- c(T.pvalue, phyper(x-1, M, N-M, n, lower.tail = F))
    k <- k+1
  }
  if(length(T.data) == 0){
    T.data <- data.frame(MCLNum = rep(1:max(T.MCLNum)),
                         T.pvalue,
                         p.adjust(T.pvalue, method = "BH"),
                         stringsAsFactors = F
    )
  }else{
    T.data <- cbind(T.data, data.frame(T.pvalue,
                                       p.adjust(T.pvalue, method = "BH"),
                                       stringsAsFactors = F
    ))
  }
  print(length(T.sample)-j)
  j <- j+1
}
colnames(T.data)[-1] <- paste0(rep(c("pvalue:", "qvalue:"), times = length(T.sample[-1])), rep(T.sample[-1], each = 2))
#save data---------------------------------------------
saveRDS(T.data, "~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_EnrichemntOfGenesList.rds")
write.table(T.data, file = "~/bigdata/yasue/PCCOfPCC/CY/Table/CY15_EnrichemntOfGenesList.txt", sep = "\t", quote = F, row.names = F)
#elapsed time------------------------------------------
after <- proc.time()
print(after - before)#104.666 sec
#remove object-----------------------------------------
rm(list = ls())