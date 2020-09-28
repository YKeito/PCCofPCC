#"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/Fig/PoP_ScatterPlot.R"
before <- proc.time()
#package-------------------------------------------------------------------------
library(dplyr)
library(stringr)
library(ggplot2)
#input data----------------------------------------------------------------------
allPCCofPCCTable <- readRDS("~/bigdata/yasue/PCCOfPCC/CY/Table/PCCofPCC/naomit/allPCCofPCCTable.rds")
T.Edge <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_MCLNumEdegeTable.rds")
T.Edge$PCC <- as.numeric(T.Edge$PCC)
T.PCC <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/PoP_CYPCC.rds")
VSPCC <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/VSPCC.rds")
#processing data-----------------------------------------------------------------
T.sample <- list.files("~/bigdata/yasue/PCCOfPCC/RDS/PublishedDataPCC/")
title <- paste0("~/bigdata/yasue/PCCOfPCC/RDS/PublishedDataPCC/", T.sample)
T.sample <- str_split(T.sample, pattern = "_PCC", simplify = T)[, 1]
T.MCLNum <- rownames(allPCCofPCCTable)[allPCCofPCCTable$Abrassicicola_qvalue < 0.05 | allPCCofPCCTable$Bcinerea_qvalue < 0.05 | 
                                         allPCCofPCCTable$Eoronti_qvalue < 0.05 | allPCCofPCCTable$MeJA_qvalue < 0.05 | 
                                         allPCCofPCCTable$PstDC3000_hrp_qvalue < 0.05 | allPCCofPCCTable$PstDC3000_qvalue < 0.05 | 
                                         allPCCofPCCTable$SA_qvalue < 0.05]
T.MCLNum <- T.MCLNum[!is.na(T.MCLNum)]
T.MCLNum <- T.MCLNum[as.numeric(T.MCLNum) < 115]
check <- list()
i <- 1
for(i in i:length(title)){
  j <- 1
  for(j in j:length(T.MCLNum)){
    T.VSPCC <- VSPCC[[i]][match(names(T.PCC[[j]]), names(VSPCC[[i]]))]
    T.which <- which(is.na(T.VSPCC))
    T.data <- allPCCofPCCTable %>% select(starts_with(T.sample[i]))
    if(length(T.which) != 0){
      check <- c(check, list(length(unique(c(str_sub(names(T.PCC[[j]][is.na(T.VSPCC)]), start = 1, end = 9),
                                             str_sub(names(T.PCC[[j]][is.na(T.VSPCC)]), start = 10, end = 18))))
                             /length(unique(c(str_sub(names(T.PCC[[j]]), start = 1, end = 9),
                                                      str_sub(names(T.PCC[[j]]), start = 10, end = 18))))
                             )
                 )
      PoP <- data.frame(CY15 = T.PCC[[j]][-T.which],
                        OpenData = T.VSPCC[-T.which]
      )
    }else{
      PoP <- data.frame(CY15 = T.PCC[[j]],
                        OpenData = T.VSPCC
      )
    }
    g <- ggplot(PoP, aes(x=CY15, y=OpenData))
    g <- g + geom_point()
    g <- g + stat_density2d(aes(fill = ..level..), alpha = 0.3, geom = "polygon")
    g <- g + scale_fill_continuous(low = "grey", high = "red", space = "Lab", name = "g = 0")
    g <- g + xlim(-1, 1) + ylim(-1, 1)
    g <- g + theme_classic()
    g <- g + xlab("CY15 PCC") + ylab(paste0(T.sample[i], " PCC"))
    g <- g + ggtitle(paste0("CY15 MCLNum", T.MCLNum[j], " : PCC = ", format(T.data[j, 1], digits = 2), ", FDR = ", format(T.data[j, 3], digits = 2)))
    g <- g + theme_bw(base_size=15)
    g <- g + theme(plot.title = element_text(hjust=0.5))
    g <- g + theme(plot.title = element_text(face = "bold"))
    g <- g + theme(axis.text=element_text(size=20))
    g <- g + theme(axis.title=element_text(size=20))
    ImageTite <- paste0("~/bigdata/yasue/PCCOfPCC/CY/Image/PoP_ScatterPlot/", T.sample[i], "/CY15MCLNum", T.MCLNum[j], "CY15PCCvs", T.sample[i], "PCC.png")
    ggsave(filename = ImageTite, plot = g, width = 7, height = 7)
    print(paste0(i, "_", length(T.MCLNum)-j))
    j <- j+1
  }
  i <- i+1
}
#save----
#saveRDS(object = T.PCC, file = "~/bigdata/yasue/PCCOfPCC/RDS/PoP_CYPCC.rds")
saveRDS(object = check, file = "~/bigdata/yasue/PCCOfPCC/RDS/PoP_remove_PCC.rds")
#elapsed time------------------------------------------
after <- proc.time()
print(after - before)#2634.408 sec 43.9068 min
#remove object-----------------------------------------------------------
rm(list = ls())