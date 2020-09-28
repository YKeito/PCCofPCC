"/home/yasue/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/DAPSeq2/UMAPCluster13_741genes_TSS.R"
#pacckage----
library(stringr)
library(dplyr)
#processing data----
before <- proc.time()
#TAIR10.gtfファイルを基にTSSとDAPSeqのピークの差を求める
TAIR10 <- read.table(file = "/root/db/TAIR10/TAIR10.gtf", fill = TRUE)
TAIR10 <- TAIR10 %>% filter(V3 == "gene")
#三つのベン図が重複する741遺伝子が対象
Target.Venn <- readRDS(file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/UMAPCluster13/UMAPCluster13_DAPSeqTargetVennTable.rds")
T.AGI <- Target.Venn %>% filter(Sample == "HY5:zfGRF:WRKYs") %>% select(AGI) %>% unlist(use.names = F)
TAIR10 <- TAIR10[match(T.AGI, TAIR10$V10), ]
#DAPSeqの結合ピーク領域のデータを獲得
DAPSeq.TFName <- readRDS("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/DAPSeq.TFName.rds")
path <- "/home/yasue/bigdata/yasue/DAPSeq/dap_data_v4/peaks/"
T.sample <- c("WRKY", "AT3G42860", "HY5")
diff.region <- c()
Direction <- c("+", "-")
chrom <- c(1:5)
f <- 1
for(f in f:length(T.sample)){
  T.Category <- DAPSeq.TFName$Category[grep(T.sample[f], DAPSeq.TFName$DAPName)]
  g <- 1
  for(g in g:length(T.Category)){
    title <- paste0(path, T.Category[g])
    names(title) <- str_split(T.Category[g], pattern = "/|_col", simplify = T)[, 2]
    Peak.data <- read.table(file = title)
    T.diff.region <- c()
    h <- 1
    for(h in h:length(chrom)){
      T.data <- Peak.data %>% filter(V1 == paste0("chr", chrom[h]))
      peak.value <- str_split(T.data$V4, pattern = ":", simplify = T)[, 2] %>% as.numeric()
      TAIR10.data <- TAIR10 %>% filter(V1 == chrom[h])
      i <- 1
      for(i in i:length(Direction)){
        direction.data <- TAIR10.data %>% filter(V7 == Direction[i])
        if(i == 1){
          TAIR10.region <- rep(direction.data$V4, each = length(peak.value))
        }else{
          TAIR10.region <- rep(direction.data$V5, each = length(peak.value))
        }
        TF.region <- rep(peak.value, times = nrow(direction.data))
        T.diff.region <- rbind(T.diff.region, 
                               data.frame(chrom = paste0("chrom", h),
                                          value = TAIR10.region-TF.region,
                                          TAIR10.region = TAIR10.region,
                                          TF.region = TF.region,
                                          direction = Direction[i],
                                          TFName = names(title),
                                          AGI = rep(direction.data$V10, each = length(peak.value)),
                                          stringsAsFactors = F)
        )
        i <- i+1
      }
      print(paste0(names(title), ":chrom", chrom[h]))
      h <- h+1
    }
    diff.region <- rbind(diff.region, T.diff.region)
    g <- g+1
  }
  f <- f+1
}
#ターゲット遺伝子に最も近いピークのみを残す
T.sample <- diff.region$TFName %>% unique()
T.AGI <- diff.region$AGI %>% unique()
summary.data <- c()
i <- 1
for(i in i:length(T.sample)){
  sample.data <- diff.region %>% filter(TFName == T.sample[i])
  j <- 1
  for(j in j:length(T.AGI)){
    T.data <- sample.data %>% filter(AGI == T.AGI[j])
    summary.data <- rbind(summary.data, T.data[which.min(abs(T.data$value)), ])
    print(paste0(length(T.sample)-i, ":", length(T.AGI)-j))
    j <- j+1
  }
  i <- i+1
}
after <- proc.time()-before
print(after)
saveRDS(object = summary.data, file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/UMAPCluster13/UMAPCluster13_741genes_TSS.rds")
rm(list = ls())