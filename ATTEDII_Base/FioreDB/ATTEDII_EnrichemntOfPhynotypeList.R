"/home/yasue/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/FioreDB/ATTEDII_EnrichemntOfPhynotypeList.R"
before <- proc.time()
#input data----
PhynotypeList <- readRDS(file = "/home/yasue/bigdata/yasue/FioreDB/PhynotypeList.rds")
MasterTable <- readRDS("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_MCLNumNodeTable.rds")
#Enrichment of Bulk Phynotype----
T.MCLNum <- max(as.numeric(MasterTable$MCLNum), na.rm = T)
allgenes <- MasterTable$AGI
N <- length(allgenes)
T.data <- c()
j <- 1
for(j in j:length(PhynotypeList)){
  T.target <- PhynotypeList[[j]]
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
                                     AGI = rep(names(PhynotypeList)[j], times = T.MCLNum),
                                     pvalue = T.pvalue,
                                     qvalue = p.adjust(T.pvalue, method = "BH"),
                                     stringsAsFactors = F
  ))
  print(length(PhynotypeList)-j)
  j <- j+1
}
#save data----
saveRDS(object = T.data, file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_EnrichemntOfPhynotypeList.rds")
#elaped time----
after <- proc.time()
print(after - before)#319.332 sec