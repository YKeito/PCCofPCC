T.sample <- list.files("~/bigdata/yasue/PCCOfPCC/RDS/PublishedDataPCC/")
title <- paste0("~/bigdata/yasue/PCCOfPCC/RDS/PublishedDataPCC/", T.sample)
T.sample <- str_split(T.sample, pattern = "_PCC", simplify = T)[, 1]
VSPCC <- list()
i <- 1
for(i in i:length(title)){
  T.VSPCCData <- readRDS(title[i])
  VSPCC <- c(VSPCC, list(T.VSPCCData$interaction_value))
  names(VSPCC[[i]]) <- paste0(T.VSPCCData$source_genes, T.VSPCCData$target_genes)
}
names(VSPCC) <- T.sample
rm(T.VSPCCData)
saveRDS(object = VSPCC, file = "~/bigdata/yasue/PCCOfPCC/RDS/VSPCC.rds")
