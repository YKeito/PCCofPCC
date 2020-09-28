"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/ATTEEDII_DAPSeqAGI.R"
pwd.family <- list.files("~/bigdata/yasue/DAPSeq/dap_data_v4/genes/")
DAPSeq_TF.AGI <- c()
i <- 1
for(i in i:length(pwd.family)){
  pwd.TF <- list.files(paste0("~/bigdata/yasue/DAPSeq/dap_data_v4/genes/", pwd.family[i], "/"))
  T.DAPSeq_TF.AGI <- c()
  j <- 1
  for(j in j:length(pwd.TF)){
    pwd <- paste0("~/bigdata/yasue/DAPSeq/dap_data_v4/genes/", pwd.family[i], "/", pwd.TF[j], "/chr1-5/chr1-5_GEM_events.nS_targets.txt")
    T.data <- read.table(pwd, , sep = "\t", quote = "", fill = T, stringsAsFactors = F, header = T)
    T.DAPSeq_TF.AGI <- c(T.DAPSeq_TF.AGI, unique(T.data$tf.at_id))
    j <- j+1
  }
  DAPSeq_TF.AGI <- rbind(DAPSeq_TF.AGI, data.frame(AGI = T.DAPSeq_TF.AGI, 
                                                   TFName = str_split(pwd.TF, pattern = "_col", simplify = T)[, 1],
                                                   pwd = paste0(rep(pwd.family[i], timess = length(T.DAPSeq_TF.AGI)), "/", pwd.TF, "/chr1-5/chr1-5_GEM_events.narrowPeak"),
                                                   stringsAsFactors = F)
                         )
  i <- i+1
}
DAPSeq_TF.AGI <- DAPSeq_TF.AGI[which(!duplicated(DAPSeq_TF.AGI$AGI)), ]
rownames(DAPSeq_TF.AGI) <- NULL
saveRDS(DAPSeq_TF.AGI, "~/bigdata/yasue/DAPSeq/DAPSeq_TF.AGI.rds")