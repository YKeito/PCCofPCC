#/home/yasue/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/DAPSeq2/UMAPCluster13_TSS_heatmap.R
#階層的クラスタリングをTSSの距離でやったバージョン
#load library----
library(dplyr)
library(tidyr)
library(ggplot2)
#WRKYsのデータを整える-----
diff.TSS <- readRDS(file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/UMAPCluster13/UMAPCluster13_741genes_TSS.rds")
WRKY.data <- diff.TSS %>% filter(grepl("WRKY", TFName))
T.AGI <- diff.TSS$AGI %>% unique()
summary.data <- c()
j <- 1
for(j in j:length(T.AGI)){
  T.data <- WRKY.data %>% filter(AGI == T.AGI[j])
  summary.data <- rbind(summary.data, T.data[which.min(abs(T.data$value)), ])
  print(length(T.AGI)-j)
  j <- j+1
}
Master.data <- rbind(diff.TSS %>% filter(TFName == "AT3G42860" | TFName == "HY5"), summary.data)
#matrix data-----
T.AGI <- Master.data$AGI %>% unique()
T.seq <- seq(-1000, 1000, 100)
T.range <- c()
i <- 1
for(i in i:c(length(T.seq)-1)){
  T.range <- c(T.range, paste0(T.seq[i], ":", T.seq[i+1]))
  i <- i+1
}
#heatmap----
T.matrix <- matrix(rep(0, times = length(T.AGI)*c(length(T.seq)-1)), nrow = length(T.AGI), ncol = c(length(T.seq)-1))
rownames(T.matrix) <- T.AGI
colnames(T.matrix) <- T.range
g <- 1
for(g in g:length(T.AGI)){
  T.data <- Master.data %>% filter(AGI == T.AGI[g])
  h <- 1
  for(h in h:c(length(T.seq)-1)){
    zfGRF <- T.seq[h] <= T.data$value[1] & T.data$value[1] < T.seq[h+1]
    if(zfGRF){
      T.matrix[g, h] <- T.matrix[g, h] + 1
    }
    HY5 <- T.seq[h] <= T.data$value[2] & T.data$value[2] < T.seq[h+1]
    if(HY5){
      T.matrix[g, h] <- T.matrix[g, h] + 2
    }
    WRKY <- T.seq[h] <= T.data$value[3] & T.data$value[3] < T.seq[h+1]
    if(WRKY){
      T.matrix[g, h] <- T.matrix[g, h] + 5
    }
    h <- h+1
  }
  print(length(T.AGI)-g)
  g <- g+1
}
CY.data <- read.table("~/Nakano_RNAseq/network_analysis/base/CY_base/allRNASeq_union.txt", sep = "\t", header = T)
CY.DEGs <- rep("No", times = nrow(T.matrix))
temp <- match(rownames(CY.data), rownames(T.matrix))
temp <- temp[!is.na(temp)]
CY.DEGs[temp] <- "Yes"
T.data <- data.frame(AGI = rep(rownames(T.matrix), times = ncol(T.matrix)),
                     value = T.matrix %>% as.vector(),
                     CY.DEGs = rep(CY.DEGs, times = ncol(T.matrix)),
                     range = rep(paste0(formatC(1:length(T.range), width = 2, flag = "0"), T.range), each = nrow(T.matrix)),
                     stringsAsFactors = F,
                     row.names = NULL)
#T.data$value[T.data$value != 3 & T.data$value != 6 & T.data$value != 8] <- 0
#T.matrix[T.matrix != 3 & T.matrix != 6 & T.matrix != 8] <- 0
T.data$value <- as.factor(T.data$value)
res <- hclust(dist(T.matrix), method = "ward.D2")
T.data$AGI <- factor(x = T.data$AGI,
                     levels = res$labels[res$order],
                     ordered = TRUE)
#colorpanel <- c("#FFFFFF", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colorpanel <- c("#FFFFFF", "#FF6666", "#56B4E9", "#CC79A7", "#009E73", "#993300", "#00FFFF", "#F0E442")
g <- ggplot(T.data, aes(x = range, y = AGI, fill = value))
g <- g + geom_tile(color = "black", size = 0.01)
g <- g + theme_bw()
#g <- g + scale_fill_gradient2(high="red", low="white")
g <- g + facet_grid(CY.DEGs~., scales="free_y", space="free", switch = "y")
g <- g + scale_fill_manual(values = colorpanel)
plot(g)
ggsave(filename = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/Image/UMAPCluster13_Fig/CommonTarget/TSS_741genes_zfGRFHeatmap.png", plot = g, width = 8, height = 6)
#heatmap_Exp
CY.data <- readRDS("~/Nakano_RNAseq/network_analysis/.RData/RDS/AverageOfallCYFPKM.rds")
CY.data <- CY.data %>% as.data.frame()
allfoldchange <- data.frame(AGI = rownames(CY.data), row.names = NULL)
T.time <- c("1h", "3h", "12h", "24h")
i <- 1
for(i in i:length(T.time)){
  DMSO.data <- CY.data %>% select(contains("DMSO")) %>% select(ends_with(T.time[i])) %>% +1 %>% as.matrix()
  CY.Exp <- CY.data %>% select(contains("CY")) %>% select(ends_with(T.time[i])) %>% +1
  allfoldchange <- cbind(allfoldchange,
                         CY.Exp %>% sweep(STATS = DMSO.data, MARGIN = 1, FUN = "/") %>% log2()
  )
  i <- i+1
}
colnames(allfoldchange)[2:ncol(allfoldchange)] <- paste0(rep(c("CY15", "CY16", "CY20"), each = 4), "_", rep(paste0(formatC(c(1, 3, 12, 24), width = 2, flag = "0"), "h"), times = 3))
T.DEGs <- rep("No", times = nrow(T.matrix)*12)
T.allfoldchange <- allfoldchange[match(rownames(T.matrix), allfoldchange$AGI), ] %>% gather(key = AGI)
colnames(T.allfoldchange)[1] <- "sample"
T.data <- T.allfoldchange %>% mutate(AGI = rep(rownames(T.matrix), times = 12),
                                     CY.Group = str_split(T.allfoldchange$sample, pattern = "_", simplify = T)[, 1],
                                     Time.Group = str_split(T.allfoldchange$sample, pattern = "_", simplify = T)[, 2],
                                     DEGs = T.DEGs
)
CY.data <- read.table("~/Nakano_RNAseq/network_analysis/base/CY_base/allRNASeq_union.txt", sep = "\t", header = T)
T.AGI <- intersect(rownames(T.matrix), rownames(CY.data))
i <- 1
for(i in i:length(T.AGI)){
  T.data$DEGs[T.data$AGI == T.AGI[i]] <- "Yes"
}
T.data$AGI <- factor(x = T.data$AGI,
                     levels = res$labels[res$order],
                     ordered = TRUE)
g <- ggplot(T.data, aes(x = sample, y = AGI, fill = value))
g <- g + geom_tile(color = "black", size = 0.01)
g <- g + theme_bw()
g <- g + scale_fill_gradient2(limits = c(-2,2), low = "blue", high = "red", na.value = "white")
g <- g + facet_grid(DEGs~CY.Group, scales="free", space="free", switch = "y")
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot(g)
ggsave(filename = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/Image/UMAPCluster13_Fig/CommonTarget/TSS_741GenesExpression_zfGRFHeatmap.png", plot = g, width = 8, height = 6)