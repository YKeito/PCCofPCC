allPCCofPCCTable <- readRDS("~/bigdata/yasue/PCCOfPCC/CY/Table/PCCofPCC/naomit/allPCCofPCCTable.rds")
T.Edge <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_MCLNumEdegeTable.rds")
#processing data-----------------------------------------------------------------
T.MCLNum <- rownames(allPCCofPCCTable)[allPCCofPCCTable$Abrassicicola_qvalue < 0.05 | allPCCofPCCTable$Bcinerea_qvalue < 0.05 | 
                                         allPCCofPCCTable$Eoronti_qvalue < 0.05 | allPCCofPCCTable$MeJA_qvalue < 0.05 | 
                                         allPCCofPCCTable$PstDC3000_hrp_qvalue < 0.05 | allPCCofPCCTable$PstDC3000_qvalue < 0.05 | 
                                         allPCCofPCCTable$SA_qvalue < 0.05]
T.MCLNum <- T.MCLNum[!is.na(T.MCLNum)]
T.MCLNum <- T.MCLNum[as.numeric(T.MCLNum) < 115]
T.PCC <- list()
j <- 1
for(j in j:length(T.MCLNum)){
  T.data <- T.Edge %>% filter(MCLNum == T.MCLNum[j])
  T.PCC <- c(T.PCC, list(T.data$PCC))
  names(T.PCC[[j]]) <- paste0(T.data$Source, T.data$Target)
}

saveRDS(object = T.PCC, file = "~/bigdata/yasue/PCCOfPCC/RDS/PoP_CYPCC.rds")