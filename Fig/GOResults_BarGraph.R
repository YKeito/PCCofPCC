"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/Fig/GOResults_BarGraph.R"
before <- proc.time()
#package----
library(ggplot2)
#input data----
allPCCofPCCTable <- readRDS("~/bigdata/yasue/PCCOfPCC/CY/Table/PCCofPCC/naomit/allPCCofPCCTable.rds")
allGOTable <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_allGOTable.rds")
#processing data----
T.MCLNum <- rownames(allPCCofPCCTable)
T.MCLNum <- T.MCLNum[as.numeric(T.MCLNum) < 115]
T.MCLNum <- intersect(T.MCLNum, allGOTable$MCLNum)
j <- 1
for(j in j:length(T.MCLNum)){
  T.data <- allGOTable[allGOTable$MCLNum == T.MCLNum[j], ]
  T.data <- T.data[order(T.data$qvalue, decreasing = F), ]
  T.data <- T.data[1:3, ]
  T.data$GOTerm <- paste0(formatC(3:1, width = 2, flag = "0"), T.data$GOTerm)
  T.data$qvalue <- -log2(T.data$qvalue)
  T.data[is.na(T.data)] <- 0
  g <- ggplot(T.data, aes(x = GOTerm, y = qvalue))
  g <- g + theme_bw()
  g <- g + geom_bar(stat = "identity")
  g <- g + coord_flip()
  g <- g + ylab("-log2(q-value)")
  g <- g + xlab("GO term")
  g <- g + theme(axis.text.x = element_text(size = 14))
  g <- g + theme(axis.text.y = element_text(size = 8))
  g <- g + theme(axis.title.x = element_blank())
  g <- g + theme(axis.title.y = element_blank())
  print(g)
  ggsave(file = paste0("~/bigdata/yasue/PCCOfPCC/CY/Image/GOResults_BarGraph/", "CY15MCLNum", T.MCLNum[j], "GOResults.png"), plot = g, width = 6, height = 2.5)
  j <- j+1
}
#elapsed time----
after <- proc.time()
print(after - before)#2.282 sec
#remove object----
rm(list = ls())