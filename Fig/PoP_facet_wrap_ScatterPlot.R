#"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/Fig/PoP_facet_wrap_ScatterPlot.R"
before <- proc.time()
#package-------------------------------------------------------------------------
library(dplyr)
library(stringr)
library(ggplot2)
#input data----------------------------------------------------------------------
allPCCofPCCTable <- readRDS("~/bigdata/yasue/PCCOfPCC/CY/Table/PCCofPCC/naomit/allPCCofPCCTable.rds")
T.Edge <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_MCLNumEdegeTable.rds")
T.PCC <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/PoP_CYPCC.rds")
VSPCC <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/VSPCC.rds")
#processing data-----------------------------------------------------------------
T.Edge$PCC <- as.numeric(T.Edge$PCC)
T.sample <- list.files("~/bigdata/yasue/PCCOfPCC/RDS/PublishedDataPCC/")
title <- paste0("~/bigdata/yasue/PCCOfPCC/RDS/PublishedDataPCC/", T.sample)
T.sample <- str_split(T.sample, pattern = "_PCC", simplify = T)[, 1]
T.MCLNum <- rownames(allPCCofPCCTable)[allPCCofPCCTable$Abrassicicola_qvalue < 0.05 | allPCCofPCCTable$Bcinerea_qvalue < 0.05 | 
                                         allPCCofPCCTable$Eoronti_qvalue < 0.05 | allPCCofPCCTable$MeJA_qvalue < 0.05 | 
                                         allPCCofPCCTable$PstDC3000_hrp_qvalue < 0.05 | allPCCofPCCTable$PstDC3000_qvalue < 0.05 | 
                                         allPCCofPCCTable$SA_qvalue < 0.05]
T.MCLNum <- T.MCLNum[!is.na(T.MCLNum)]
T.MCLNum <- T.MCLNum[as.numeric(T.MCLNum) < 115]
j <- 1
for(j in j:length(T.MCLNum)){
  PoP <- c()
  i <- 1
  for(i in i:length(title)){
    T.VSPCC <- VSPCC[[i]][match(names(T.PCC[[j]]), names(VSPCC[[i]]))]
    T.which <- which(is.na(T.VSPCC))
    T.data <- allPCCofPCCTable %>% select(starts_with(T.sample[i]))
    if(length(T.which) != 0){
      PoP <- rbind(PoP, 
                   data.frame(CY15 = T.PCC[[j]][-T.which],
                              OpenData = T.VSPCC[-T.which],
                              sample = rep(paste0(T.sample[i], " : PCC = ", format(T.data[j, 1], digits = 2), ", FDR = ", format(T.data[j, 3], digits = 2)),
                                           times = length(T.VSPCC[-T.which]))
                   )
      )
    }else{
      PoP <- rbind(PoP, 
                   data.frame(CY15 = T.PCC[[j]],
                              OpenData = T.VSPCC,
                              sample = rep(paste0(T.sample[i], ":PCC=", format(T.data[j, 1], digits = 2), ", FDR=", format(T.data[j, 3], digits = 2)), 
                                           times = length(T.VSPCC)))
      )
    }
    i <- i+1
  }
  g <- ggplot(PoP, aes(x=CY15, y=OpenData))
  g <- g + geom_point()
  g <- g + xlab("CY15 PCC") + ylab("OpenData PCC")
  g <- g + facet_wrap(~ sample, ncol = 4)
  g <- g + xlim(-1, 1) + ylim(-1, 1)
  g <- g + theme_classic()
  #g <- g + stat_density2d(aes(alpha=..density..), geom = "tile",contour = FALSE)
  g <- g + stat_density2d(aes(fill = ..level..), alpha = 0.3, geom = "polygon")
  g <- g + scale_fill_continuous(low = "grey", high = "red", space = "Lab", name = "g = 0")
  g <- g +  scale_colour_discrete(guide = FALSE)
  g <- g + theme(strip.text = element_text(size = 9))
  g <- g + theme(axis.text.x = element_text(size = 10))
  g <- g + theme(axis.text.y = element_text(size = 10))
  ImageTite <- paste0("~/bigdata/yasue/PCCOfPCC/CY/Image/PoP_ScatterPlot/facet_wrap/CY15MCLNum", T.MCLNum[j], "CY15PCCvsOpenDataPCC.png")
  ggsave(filename = ImageTite, plot = g, width = 12, height = 6)
  #save
  RDStitle <- paste0("~/bigdata/yasue/PCCOfPCC/RDS/PoP_facet_wrap/CY15MCLNum", T.MCLNum[j], "CY15PCCvsOpenData.rds")
  saveRDS(object = PoP, file = RDStitle)
  print(length(T.MCLNum)-j)
  j <- j+1
}
#elapsed time------------------------------------------
after <- proc.time()
print(after - before)#4234.096 sec, 70.56827 min
#remove object-----------------------------------------------------------
rm(list = ls())