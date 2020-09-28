#"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/CY15_EdgeInfo.R"
before <- proc.time()
#package--------------------
library(stringr)
#input Data-----------------
EdgeTable <- read.table("~/bigdata/yasue/PCCOfPCC/CY/Table/NetworkAnalyzerResults/CY15EdgeTable_default.csv", sep = ",", header = T, stringsAsFactors = F)
NodeTable <- read.table("~/bigdata/yasue/PCCOfPCC/CY/Table/NetworkAnalyzerResults/CY15NodeTable_default.csv", sep = ",", header = T, stringsAsFactors = F)
#processing data------------
T.Edge <- str_split(EdgeTable$name, pattern = " ", simplify = T)
T.Source <- T.Edge[, 1]
T.Target <- T.Edge[, 3]
T.MCLNum <- NodeTable$X__mclCluster
names(T.MCLNum) <- NodeTable$name
T.data <- data.frame(Source = T.Source,
                     Target = T.Target,
                     MCLNum = rep(0, times = length(T.Target)),
                     PCC = str_sub(T.Edge[, 2], start = 2, end = c(nchar(T.Edge[, 2])-1)),
                     stringsAsFactors = F
)
n <- 1
for(n in n:length(T.MCLNum)){
  T.data$MCLNum[names(T.MCLNum)[n] == T.data$Source] <- T.MCLNum[n]
  print(length(T.MCLNum)-n)
  n <- n+1
}
#save data------------------
saveRDS(T.data, "~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_MCLNumEdegeTable.rds")
#elapsed time------------------------------------------
after <- proc.time()
print(after - before)#117.647 sec
#remove object---------------
rm(list = ls())