#"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/Motif_EnrichmentOfCY15.R"
before <- proc.time()
#package---------------------------------------------------------
library(stringr)
library(dplyr)
#input data------------------------------------------------------
MasterTable <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_MCLNumNodeTable.rds")
CY15_EnrichemntOfGenesList <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_EnrichemntOfGenesList.rds")
upstream_500 <- read.table(file = "~/bigdata/yasue/motif/TAIR10_upstream_500_20101028", fill = T, sep = ",", stringsAsFactors = F)
upstream_500 <- as.character(unlist(upstream_500))
#processing data-------------------------------------------------
temp <- grep("chr", upstream_500)
T.AGI <- c()
i <- 1
total <- length(temp)
for(i in i:total){
  T.AGI <- c(T.AGI, substr(upstream_500[temp[i]], 2, 10))
  print(i)
  i <- i+1
}
#最後の一周だけ自動化できなかった。
allsequence <- list()
presequence <- c()
sequence <- c()
i <- 1
total <- length(temp)
for(i in i:c(total-1)){
  presequence <- upstream_500[c(temp[i]+1):c(temp[i+1]-1)]
  n <- 1
  for(n in n:length(presequence)){
    sequence <- paste0(sequence, presequence[n])
    n <- n+1
  }
  allsequence <- c(allsequence, list(sequence))
  sequence <- c()
  print(i)
  i <- i+1
}
#残りの最後の一周を追加
presequence <- upstream_500[c(temp[total]+1):length(upstream_500)]
sequence <- paste0(sequence, presequence[n])
allsequence <- c(allsequence, list(sequence))
#data.frame
up_500bp <- data.frame(AGI = T.AGI,
                       sequence = unlist(allsequence),
                       stringsAsFactors = F
)
#各MCLNumのAGIの配列をクラスター単位で引っこ抜く------------------------------------------------------
T.data <- CY15_EnrichemntOfGenesList %>% select("MCLNum", starts_with("qvalue:CY15"))
T.MCLNum <- T.data[, 1][T.data[, 2] < 0.05 | T.data[, 3] < 0.05 | T.data[, 4] < 0.05 | T.data[, 5] < 0.05]
i <- 1
for(i in i:length(T.MCLNum)){
  T.AGI <- MasterTable %>% filter(MCLNum == T.MCLNum[i]) %>% select("AGI") %>% unlist()
  T.data <- up_500bp[match(T.AGI, up_500bp$AGI), ]
  T.Control <- sample(up_500bp$AGI, length(T.AGI))
  T.ControlData <- up_500bp[match(T.Control, up_500bp$AGI), ]
  
  fasta <- c()
  allfasta <- c()
  cont_fasta <- c()
  cont_allfasta <- c()
  total_o <- nrow(T.data)
  o <- 1
  for(o in o:total_o){
    fasta <- rbind(str_sub(T.data[, "sequence"][o], start=1, end=80),
                   str_sub(T.data[, "sequence"][o], start=81, end=160),
                   str_sub(T.data[, "sequence"][o], start=161, end=240),
                   str_sub(T.data[, "sequence"][o], start=241, end=320),
                   str_sub(T.data[, "sequence"][o], start=321, end=400),
                   str_sub(T.data[, "sequence"][o], start=401, end=480),
                   str_sub(T.data[, "sequence"][o], start=481, end=500)
    )
    cont_fasta <- rbind(str_sub(T.ControlData[, "sequence"][o], start=1, end=80),
                        str_sub(T.ControlData[, "sequence"][o], start=81, end=160),
                        str_sub(T.ControlData[, "sequence"][o], start=161, end=240),
                        str_sub(T.ControlData[, "sequence"][o], start=241, end=320),
                        str_sub(T.ControlData[, "sequence"][o], start=321, end=400),
                        str_sub(T.ControlData[, "sequence"][o], start=401, end=480),
                        str_sub(T.ControlData[, "sequence"][o], start=481, end=500)
    )
    data_fastaAGI <- paste0(">", T.data[, "AGI"][o])
    control_fastaAGI <- paste0(">", T.ControlData[, "AGI"][o])
    allfasta <- c(allfasta, rbind(data_fastaAGI, fasta))
    cont_allfasta <- c(cont_allfasta, rbind(control_fastaAGI, cont_fasta))
    o <- o+1
  }
  target <- paste0("~/bigdata/yasue/PCCOfPCC/Motif/MultiFasta/Target/", "MCLNum", T.MCLNum[i], "_upstream500.fasta")
  control <- paste0("~/bigdata/yasue/PCCOfPCC/Motif/MultiFasta/Control/", "control_MCLNum", T.MCLNum[i], "_upstream500.fasta")
  write.table(allfasta, file = target, append = F, quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(cont_allfasta, file = control, append = F, quote = F, sep = "\t", row.names = F, col.names = F)
  print(i)
  i <- i+1
}
#elapsed time-----------------------------------------------------
after <- proc.time()
print(after - before)#61.758 sec
#remove object----------------------------------------------------
rm(list = ls())