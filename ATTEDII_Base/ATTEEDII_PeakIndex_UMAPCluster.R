"/home/yasue/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/DAPSeq/ATTEEDII_PeakIndex_UMAPCluster.R"
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
    TFName <- unique(as.vector(temp))
    T.matrix <- matrix(rep(1, times = length(TFName)*length(TFName)), nrow = length(TFName), ncol = length(TFName))
    rownames(T.matrix) <- TFName
    colnames(T.matrix) <- TFName
    i <- 1
    for(i in i:length(T.filename)){
      pwd <- paste0(path, UMAPDir[h], "/", T.filename[i])
      T.data <- readRDS(file = pwd)
      T.value <- T.data %>% filter(value <= 100 & value >= -100) %>% nrow()/T.data %>% nrow()
      T.matrix[which(rownames(T.matrix) == temp[i, 1]), which(colnames(T.matrix) == temp[i, 2])] <- T.value
      T.matrix[which(rownames(T.matrix) == temp[i, 2]), which(colnames(T.matrix) == temp[i, 1])] <- T.value
      #print(paste0(UMAPDir[h], ":", length(T.filename)-i))
      i <- i+1
    }
    res <- hclust(dist(T.matrix), method = "ward.D2")
    df <- T.matrix %>% melt(id.var="AF") %>% mutate(Cluster = rep(UMAPDir[h], times = length(T.value)), diffTFfamily = rep(0, times = length(T.value)))
    T.names <- permutations(n = length(TFName), r = 2, v=TFName, repeats.allowed = T)
    diffTFfamily <- c()
    i <- 1
    for(i in i:nrow(T.names)){
      diffTFfamily <- rbind(diffTFfamily, data.frame(Var1 = paste0(T.names[i, 1], T.names[i, 2]), 
                                                     diff = T.family[T.names[i, 1] == names(T.family)] == T.family[T.names[i, 2] == names(T.family)],
                                                     row.names = NULL)
                            )
      i <- i+1
    }
    df$diffTFfamily[match(diffTFfamily$Var1, paste0(df$Var1, df$Var2))] <- as.character(diffTFfamily$diff)
    #library("gplots")
    #heatmap.2(T.matrix, scale = "none", col = bluered(100), trace = "none", density.info = "none", breaks = seq(0, 1, 0.01))
    df$Var1 <- factor(x = df$Var1,
                            levels = res$labels[res$order],
                            ordered = TRUE)
    df$Var2 <- factor(x = df$Var2,
                      levels = res$labels[res$order],
                      ordered = TRUE)
    PeakIndex <- rbind(PeakIndex, df)
    heatmap.plot <- ggplot(NULL)
    heatmap.plot <- heatmap.plot + geom_tile(data = df, aes(x = Var1, y = Var2, fill = value), color = "white", size = 0.05)
    heatmap.plot <- heatmap.plot + geom_tile(data = df %>% filter(value > 0.5 & diffTFfamily == FALSE), aes(x = Var1, y = Var2, fill = value), color = "black", size = 0.2)
    heatmap.plot <- heatmap.plot + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
    heatmap.plot <- heatmap.plot + theme(axis.text = element_text(size = 8))
    heatmap.plot + theme(legend.position = 'none')
    #heatmap.plot <- heatmap.plot + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
    heatmap.plot <- heatmap.plot + scale_fill_gradient(low = "white", high = "red", name="Volume", limits = c(0, 1))
    heatmap.plot <- heatmap.plot + ggtitle(paste0(UMAPDir[h], ":PeakIndex_Heatmap"))
    heatmap.plot <- heatmap.plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    heatmap.plot <- heatmap.plot + theme(title = element_text(size = 5))
    heatmap.plot <- heatmap.plot + theme(legend.position = 'none')
    plot(heatmap.plot)
    ggsave(paste0("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/Image/20190809DAPSeq_GenomeWide/PeakIndex/100bp/PeakIndex100bp_Heatmap", UMAPDir[h], ".png"), heatmap.plot, width = 10, height = 5)
  }
  h <- h+1
}
saveRDS(object = PeakIndex, file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/PeakIndex100bp.rds")
#check
PeakIndex %>% filter(value > 0.5 & diffTFfamily == FALSE)#0.15とかまでありかも0.1はびみょい
PeakIndex %>% filter(Var1 == "MYB3R1_col_a" | Var1 == "MYB3R4_col_a" | Var1 == "MYB3R5_col_a") %>% filter(Var2 == "MYB3R1_col_a" | Var2 == "MYB3R4_col_a" | Var2 == "MYB3R5_col_a")
#elapsed time----
after <- proc.time()
print(after - before)#167.525 sec