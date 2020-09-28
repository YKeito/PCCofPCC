#"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/CY15MCLNum_EnrichedMotif.R"
before <- proc.time()
#input data------------------------------------------------
T.filename <- list.files("~/bigdata/yasue/PCCOfPCC/Motif/MultiFasta/Results", full.names = T)
CY15_EnrichemntOfGenesList <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_EnrichemntOfGenesList.rds")
#processing data----------
T.data <- CY15_EnrichemntOfGenesList[CY15_EnrichemntOfGenesList$MCLNum < 115, ]
T.data <- T.data[T.data$`qvalue:CY15_12h_FDR0.05` < 0.05 | T.data$`qvalue:CY15_1h_FDR0.05` < 0.05 |
                   T.data$`qvalue:CY15_24h_FDR0.05` < 0.05 | T.data$`qvalue:CY15_3h_FDR0.05` < 0.05 |
                   T.data$`qvalue:Botrytis_cinerea` < 0.05 | T.data$`qvalue:PstDC3000_WT` < 0.05 | 
                   T.data$`qvalue:PstDC3000_hrp` < 0.05 | T.data$`qvalue:MeJA_DEGs` < 0.05 | 
                   T.data$`qvalue:SA_DEGs` < 0.05, ]
T.data <- T.data %>% select(-starts_with("pvalue"))
T.data <- T.data %>% select(c(-contains("BTH"), -contains("CY16"), -contains("CY20"), -contains("MCLNum")))
T.MCLNum <- rownames(T.data)
MotifID <- c()
obnames <- c()
consensus <- c()
motif_alt_ID <- c()
i <- 1
for(i in i:length(T.MCLNum)){
  filename <- T.filename[grep(paste0("MCLNum", T.MCLNum[i], "_"), T.filename)]
  title <- paste0(filename, "/ame.tsv")
  motif_ame_results <- try(read.table(title, sep = "\t", header = T, stringsAsFactors = F), silent = TRUE)
  if(class(motif_ame_results) != "try-error"){
    MotifID <- c(MotifID, motif_ame_results$motif_ID)
    consensus <- c(consensus, motif_ame_results$consensus)
    motif_alt_ID <- c(motif_alt_ID, motif_ame_results$motif_alt_ID)
    obnames <- c(obnames, rep(T.MCLNum[i], times = length(motif_ame_results$motif_ID)))
  }
  i <- i+1
}
MotifID_consensus <- data.frame(MCLNum = obnames,
                                     MotifID = MotifID,
                                     consensu = consensus,
                                     motif_alt_ID = motif_alt_ID,
                                     stringsAsFactors = F
)
#save data------
saveRDS(MotifID_consensus, "~/bigdata/yasue/PCCOfPCC/RDS/Table/CT15EnrichedMotif.rds")
write.table(MotifID_consensus, file = "~/bigdata/yasue/PCCOfPCC/CY/Table/MotifResults/CT15EnrichedMotif.txt", sep = "\t", quote = F, row.names = F)
#elapsed time------------------------------------------
after <- proc.time()
print(after - before)#2889.739 sec, 48.16232 min
#remove object-----------------------------------------
rm(list = ls())
