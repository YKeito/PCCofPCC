#"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/ATTEDII_GOHeatmap.R"
#package----
library(ggplot2)
library(stringr)
library(gplots)
library(RColorBrewer)
#input data----
allGOTable <- readRDS(file = "~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/ATTEDII_allGOTable.rds")
allGOTable <- allGOTable %>% select(-"pvalue") 
allGOTable$MCLNum <- formatC(allGOTable$MCLNum, width = 3, flag = 0)
#Heatmap----
#MCLNum
allGOTable$qvalue[allGOTable$qvalue > 5e-2] <- 400
allGOTable$qvalue[allGOTable$qvalue <= 5e-2 & allGOTable$qvalue > 5e-5] <- 300
allGOTable$qvalue[allGOTable$qvalue <= 5e-5 & allGOTable$qvalue > 5e-10] <- 200
allGOTable$qvalue[allGOTable$qvalue <= 5e-10 & allGOTable$qvalue > 5e-20] <- 100
allGOTable$qvalue[allGOTable$qvalue <= 5e-20] <- 0
#Create heatmap plot
heatmap.plot <- ggplot(data = allGOTable, aes(x = MCLNum, y = GOTerm, fill = qvalue))
heatmap.plot <- heatmap.plot +geom_tile(color = "black", size = 0.1)
heatmap.plot <- heatmap.plot + theme_classic()
heatmap.plot <- heatmap.plot + scale_fill_gradient(low="seagreen",high="white")
heatmap.plot <- heatmap.plot + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
heatmap.plot <- heatmap.plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
heatmap.plot <- heatmap.plot + theme(legend.position = 'none')
#heatmap.plot <- heatmap.plot + theme(axis.text=element_text(size=1))
#save data----
ggsave("~/bigdata/yasue/PCCOfPCC/ATTEDII/Image/EnrichmentOfGenesList/ATTEDII_GOHeatmap.png", heatmap.plot, width = 24, height = 8)
#remove object----
rm(list = ls())