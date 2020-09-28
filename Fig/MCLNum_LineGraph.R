"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/Fig/MCLNum_LineGraph.R"
before <- proc.time()
#package-------------------------------------------------------------------------
library(dplyr)
library(stringr)
library(ggplot2)
#input data-----
CY.FPKM <- readRDS("~/Nakano_RNAseq/network_analysis/.RData/RDS/allCYFPKM.rds")
allPCCofPCCTable <- readRDS("~/bigdata/yasue/PCCOfPCC/CY/Table/PCCofPCC/naomit/allPCCofPCCTable.rds")
MasterTable <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_MCLNumNodeTable.rds")
#processing data----
T.MCLNum <- rownames(allPCCofPCCTable)[allPCCofPCCTable$Abrassicicola_qvalue < 0.05 | allPCCofPCCTable$Bcinerea_qvalue < 0.05 | 
                                         allPCCofPCCTable$Eoronti_qvalue < 0.05 | allPCCofPCCTable$MeJA_qvalue < 0.05 | 
                                         allPCCofPCCTable$PstDC3000_hrp_qvalue < 0.05 | allPCCofPCCTable$PstDC3000_qvalue < 0.05 | 
                                         allPCCofPCCTable$SA_qvalue < 0.05]
T.MCLNum <- T.MCLNum[!is.na(T.MCLNum)]
T.MCLNum <- T.MCLNum[as.numeric(T.MCLNum) < 115]
T.time <- c("_1h", "_3h", "_12h", "_24h")
T.condition <- paste0(rep(c("DMSO", "CY15"), each = 4), rep(T.time, times = 2))
i <- 1
for(i in i:length(T.MCLNum)){
  T.AGI <- MasterTable %>% filter(MCLNum == T.MCLNum[i]) %>% select("AGI") %>% unlist() %>% as.vector()
  T.average <- list()
  j <- 1
  for(j in j:length(T.condition)){
    T.average <- c(T.average, list(CY.FPKM[match(T.AGI, rownames(CY.FPKM)), ] %>% select(starts_with(T.condition[j])) %>% apply(MARGIN = 1, FUN = mean)))
    names(T.average)[[j]] <- T.condition[j]
    j <- j+1
  }
  T.average <- data.frame(T.average, stringsAsFactors = F)
  T.FPKM <- list()
  j <- 1
  for(j in j:length(T.time)){
    T.data <- T.average %>% select(ends_with(T.time[j]))
    T.FPKM <- c(T.FPKM, list(log2(c(T.data[, 2]+0.01)/c(T.data[, 1]+0.01))))
    names(T.FPKM[[j]]) <- rownames(T.average)
    j <- j+1
  }
  names(T.FPKM) <- paste0("CY15", T.time, "/DMSO", T.time)
  T.FPKM <- data.frame(T.FPKM, stringsAsFactors = F)
  T.FPKM <- T.FPKM[order(T.FPKM$CY15_1h.DMSO_1h, decreasing = T), ]
  df <- data.frame(AGI = rep(rownames(T.FPKM), times = 4),
                   time = rep(c(1, 3, 12, 24), each = nrow(T.average)),
                   log2FC = unlist(T.FPKM),
                   heatmap_reorder = rep(1:nrow(T.FPKM), times = 4),
                   group = rep("", times = length(unlist(T.FPKM))),
                   stringsAsFactors = F)
  j <- 1
  PCC <- c()
  for(j in j:nrow(T.FPKM)){
    PCC <- c(PCC, cor(as.numeric(T.FPKM[1, ]), as.numeric(T.FPKM[j, ])))
  }
  df$group[c(which(PCC < 0), nrow(T.FPKM)+which(PCC < 0), nrow(T.FPKM)*2+which(PCC < 0), nrow(T.FPKM)*3+which(PCC < 0))] <- paste0("02negative ", sum(PCC < 0), " genes")
  df$group[c(which(PCC > 0), nrow(T.FPKM)+which(PCC > 0), nrow(T.FPKM)*2+which(PCC > 0), nrow(T.FPKM)*3+which(PCC > 0))] <- paste0("01possitive ", sum(PCC > 0), " genes")
  #geom_smooth
  g <- ggplot(data = df, aes(x = time, y = log2FC, colour = group))
  #g <- g + geom_point(mapping = aes(x = time, y = log2FC, colour = group))
  #g <- g + geom_line(mapping = aes(x = time, y = log2FC, group = AGI))
  g <- g + geom_smooth(method = "loess", mapping = aes(x = time, y = log2FC, colour = group), level = 0.95)
  g <- g + theme(legend.position="top")
  g <- g + ggtitle(paste0("CY15 MCLNum", T.MCLNum[i]))
  g <- g + theme(axis.text.x = element_text(size = 14))
  g <- g + theme(axis.text.y = element_text(size = 14))
  title <- paste0("~/bigdata/yasue/PCCOfPCC/CY/Image/MCLNumLineGraph/CY15/CY15MCLNum", T.MCLNum[i], "_smooth", ".png")
  ggsave(title, g, width = 7, height = 5)
  #geom_tile
  df$time <- paste0(formatC(df$time, width = 2, flag = "0"), "h")
  df$log2FC[df$log2FC < -3] <- -3
  df$log2FC[df$log2FC > 3] <- 3
  g1 <- ggplot(df, aes(x = time, y = reorder(AGI, heatmap_reorder), fill = log2FC))
  g1 <- g1 + geom_tile(color = "black", size = 0.1)
  g1 <- g1 + theme_bw()
  g1 <- g1 + scale_fill_gradient2(limits = c(-3,3), low = "blue", high = "red", na.value = "white")
  g1 <- g1 + ggtitle(paste0("CY15 MCLNum", T.MCLNum[i]))
  title <- paste0("~/bigdata/yasue/PCCOfPCC/CY/Image/MCLNumHeatMap/CY15/CY15MCLNum", T.MCLNum[i], "_heatmap", ".png")
  ggsave(title, g1, height = 7, width = 5)
  i <- i+1
}
#elapsed time------------------------------------------
after <- proc.time()
print(after - before)#8117.808 sec 135.2968 min
#remove object-----------------------------------------------------------
rm(list = ls())