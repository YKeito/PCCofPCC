#"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/PCCofPCC_CY15vsOpendata"
before <- proc.time()
#package-------------------------------------------------------------------------
library(Hmisc)
library(dplyr)
library(stringr)
library(ggplot2)
#input data----------------------------------------------------------------------
MasterTable <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_MCLNumNodeTable.rds")
T.Edge <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_MCLNumEdegeTable.rds")
T.Edge$PCC <- as.numeric(T.Edge$PCC)
CY15_EnrichemntOfGenesList <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_EnrichemntOfGenesList.rds")
#processing data-----------------------------------------------------------------
T.sample <- list.files("~/bigdata/yasue/PCCOfPCC/RDS/PublishedDataPCC/")
title <- paste0("~/bigdata/yasue/PCCOfPCC/RDS/PublishedDataPCC/", T.sample)
T.sample <- str_split(T.sample, pattern = "_PCC", simplify = T)[, 1]
T.data <- CY15_EnrichemntOfGenesList %>% select("MCLNum", starts_with("qvalue:CY15"))
T.MCLNum <- T.data[, 1][T.data[, 2] < 0.05 | T.data[, 3] < 0.05 | T.data[, 4] < 0.05 | T.data[, 5] < 0.05]
check <- list()
allPCCofPCCTable <- c()
i <- 1
for(i in i:length(title)){
  PCC <- c()
  Pvalue <- c()
  j <- 1
  for(j in j:length(T.MCLNum)){
    T.VSPCC <- VSPCC[[i]][match(names(T.PCC[[j]]), names(VSPCC[[i]]))]
    T.which <- which(is.na(T.VSPCC))
    if(length(T.which) != 0){
      T.VSPCC[[j]] <- T.VSPCC[[j]][-T.which]
      T.PCC[[j]] <- T.PCC[[j]][-T.which]
    }
    if(length(T.PCC) >= 3){
      T.PCCofPCC <- cor.test(T.PCC, T.VSPCC)
      PCC <- c(PCC, T.PCCofPCC$estimate)
      names(PCC)[j] <- T.MCLNum[j]
      Pvalue <- c(Pvalue, T.PCCofPCC$p.value)
      names(Pvalue)[j] <- T.MCLNum[j]
    }else{
      PCC <- c(PCC, NA)
      names(PCC)[j] <- T.MCLNum[j]
      Pvalue <- c(Pvalue, NA)
      names(Pvalue)[j] <- T.MCLNum[j]
    }
    print(paste0(i, "_", length(T.MCLNum)-j))
    j <- j+1
  }
  PCCofPCCTable <- data.frame(MCLNum = T.MCLNum,
                              PoP = PCC,
                              pvalue = Pvalue,
                              qvalue = p.adjust(Pvalue, method = "BH"),
                              row.names = 1,
                              stringsAsFactors = F
  )
  check <- c(check, list(rownames(PCCofPCCTable)[PCCofPCCTable$qvalue < 0.05]))
  colnames(PCCofPCCTable) <- paste0(T.sample[i], "_", colnames(PCCofPCCTable))
  if(length(allPCCofPCCTable) == 0){
    allPCCofPCCTable <- PCCofPCCTable
  }else{
    allPCCofPCCTable <- cbind(allPCCofPCCTable, PCCofPCCTable)
  }
  write.table(PCCofPCCTable, file = paste0("~/bigdata/yasue/PCCOfPCC/CY/Table/PCCofPCC/naomit/PoPTableOf", T.sample[i], "_naomit.txt"), append=F, quote = F, sep = "\t", row.names = T)
  i <- i+1
}
saveRDS(allPCCofPCCTable, "~/bigdata/yasue/PCCOfPCC/CY/Table/PCCofPCC/naomit/allPCCofPCCTable.rds")
write.table(allPCCofPCCTable, "~/bigdata/yasue/PCCOfPCC/CY/Table/PCCofPCC/naomit/allPCCofPCCTable.txt", append=F, quote = F, sep = "\t", row.names = T)
#elapsed time------------------------------------------
after <- proc.time()
print(after - before)#39868.242 sec 664.4707 min 11.07451 hr
#remove object-----------------------------------------------------------
check
rm(list = ls())