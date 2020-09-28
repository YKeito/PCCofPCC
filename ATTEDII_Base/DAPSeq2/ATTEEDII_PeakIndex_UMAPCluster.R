"/home/yasue/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/DAPSeq2/ATTEEDII_PeakIndex_UMAPCluster.R"
#pacckage----
library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(gtools)
#processing data----
before <- proc.time()
DAPSeq.TFName <- readRDS("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/DAPSeq.TFName.rds")
T.family <- str_split(DAPSeq.TFName$Category, pattern = "/", simplify = T)[, 1]
names(T.family) <- str_split(DAPSeq.TFName$Category, pattern = "/", simplify = T)[, 2]
path <- "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/"
UMAPDir <- list.files(path = path)
UMAPDir <- UMAPDir[grep("UMAPCluster", UMAPDir)]
PeakIndex <- c()
h <- 1
for(h in h:length(UMAPDir)){
  T.filename <- list.files(paste0(path, UMAPDir[h]))
  if(length(T.filename) > 1){
    temp <- str_split(str_split(T.filename, pattern = "DAPSeq_base", simplify = T)[, 2], pattern = "&", simplify = T)
    temp[, 2] <- str_split(str_split(temp[, 2], pattern = "target", simplify = T)[, 2], pattern = ".rds", simplify = T)[, 1]
    T.PeakIndex <- c()
    i <- 1
    for(i in i:length(T.filename)){
      pwd <- paste0(path, UMAPDir[h], "/", T.filename[i])
      T.data <- readRDS(file = pwd)
      hist <- hist(T.data$value, breaks = seq(-3000, 3000, 50))
      T.value <- rbind(hist$counts)
      colnames(T.value) <- paste0("PeakCount:", hist$mids)
      T.PeakIndex <- rbind(T.PeakIndex, data.frame(source = temp[i, 1],
                                                   sAGI = DAPSeq.TFName$AGI[grep(temp[i, 1], names(T.family))],
                                                   target = temp[i, 2],
                                                   tAGI = DAPSeq.TFName$AGI[grep(temp[i, 2], names(T.family))],
                                                   TotalPeaks = nrow(T.data),
                                                   Cluster = UMAPDir[h],
                                                   diffTFfamily = as.character(T.family[temp[i, 1] == names(T.family)] == T.family[temp[i, 2] == names(T.family)]),
                                                   PeakMax = max(hist$counts),
                                                   Peak_which.max = hist$mids[which.max(hist$counts)],
                                                   T.value,
                                                   row.names = NULL,
                                                   stringsAsFactors = F)
      )
      print(paste0(UMAPDir[h], ":", length(T.filename)-i))
      i <- i+1
    }
    PeakIndex <- rbind(PeakIndex, T.PeakIndex)
  }
  h <- h+1
}
Bin.max <- PeakIndex %>% select(contains("PeakCount")) %>% apply(MARGIN = 1, FUN = which.max) #最も平均Peakが高い結合領域はどこか
th_200bp <- PeakIndex %>% select("PeakCount..175":"PeakCount.175") %>% apply(MARGIN = 1, FUN = sum) %>% as.matrix() %>% sweep(MARGIN = 1, STATS = 400, FUN = "/")
th_3000bp <- PeakIndex %>% select(contains("PeakCount"), -c("PeakCount..175":"PeakCount.175")) %>% apply(MARGIN = 1, FUN = sum) %>% as.matrix() %>% sweep(MARGIN = 1, STATS = 5600, FUN = "/")
PeakIndex <- PeakIndex %>% mutate(Bin_which.max = Bin.max, PeakIndex = th_200bp/th_3000bp)
saveRDS(object = PeakIndex, file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/PeakIndex.rds")
test.data <- test.data %>% select(-contains("PeakCount"))
test.data <- test.data %>% filter(PeakIndex > 2) #4578
test.data <- test.data %>% subset(PeakMax > 50) #2578
saveRDS(object = test.data, file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/PeakIndex_thver.rds")
#chechk:histogram----
t.test.data <- PeakIndex %>% mutate(Bin_which.max = Bin.max, PeakIndex = th_200bp/th_3000bp) %>% select(-contains("PeakCount")) %>% filter(Cluster == "UMAPCluster13") %>% filter(diffTFfamily == FALSE) %>% filter(PeakIndex > 2) %>% subset(PeakMax > 50)
g <- 1
for(g in g:nrow(t.test.data)){
  title <- paste0("/home/yasue//bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/", t.test.data$Cluster[g], "/DAPSeq_base", t.test.data$source[g], "&target", t.test.data$target[g], ".rds")
  allchrom.LengthDifference <- readRDS(title)
  df <- ggplot(allchrom.LengthDifference, aes(x = value))
  df <- df + geom_histogram(position="identity", alpha = 0.6, binwidth = 50)
  df <- df + theme_bw()
  df <- df + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  df <- df + theme(legend.position = 'none')
  df <- df + theme(axis.text=element_text(size=18))
  df <- df + scale_x_continuous(breaks = seq(-3000, 3000, 1000))
  df <- df + ggtitle(paste0(t.test.data$source[g], t.test.data$target[g]))
  plot(df)
  ggsave(filename = paste0("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/Image/UMAPCluster13_Fig/", "DAPSeq_base", t.test.data$source[g], "&target", t.test.data$target[g], ".png"), 
         plot = df, width = 8, height = 6)
  g <- g+1
}
#elapsed time----
after <- proc.time()
print(after - before)#601.733 sec