#"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/allSAFPKM.R"
library(dplyr)
library(stringr)
pwd <- list.files("~/bigdata/yasue/OpenData/PRJNA224133/cufflinks", full.names = T)
cname <- str_split(pwd, pattern = "cufflinks/", simplify = T)[, 2]
T.data <- c()
n <- 1
for(n in n:length(pwd)){
  title <- paste0(pwd[n], "/genes.fpkm_tracking")
  T.FPKM <- read.table(title, header = T, stringsAsFactors = F, fill = T, sep = "\t", quote = "")
  T.FPKM <- T.FPKM %>% select("tracking_id", "FPKM")
  if(is.null(T.data)){
    T.data <- T.FPKM
  }else{
    T.data <- cbind(T.data, T.FPKM$FPKM)
  }
  print(n)
}
colnames(T.data) <- c("AGI", cname)
T.data <- T.data[!duplicated(T.data$AGI), ]
rownames(T.data) <- T.data$AGI
T.data <- T.data[, 2:ncol(T.data)]
saveRDS(T.data, "~/bigdata/yasue/GEO/RDS/allSAFPKM.rds")