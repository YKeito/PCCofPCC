"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/Fig/CY15andGenesList_heatmap.R"
#package----
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(gplots)
#input data----
CY15_EnrichemntOfGenesList <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_EnrichemntOfGenesList.rds")
Hcluster.Results <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_MCLNumDendrogram.rds")
SHcluster.Results <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_SampleDendrogram.rds")
#processing data----
T.data <- CY15_EnrichemntOfGenesList %>% select(-starts_with("pvalue"), -contains("BTH"), -contains("CY16"), -contains("CY20"))
T.data <- T.data[apply(T.data < 0.05, MARGIN = 1, FUN = sum) != 0, ]
T.data <- T.data %>% gather(key = Sample, value = FDR, -MCLNum)
#Heatmap----
T.data$MCLNum <- factor(x = T.data$MCLNum,
                        levels = Hcluster.Results$MCLNum[Hcluster.Results$order],
                        ordered = TRUE)
T.sample <- str_split(unique(T.data$Sample), pattern = ":", simplify = T)[, 2]
T.data$Sample <- rep(T.sample, each = length(unique(T.data$MCLNum)))
T.data$Sample <- factor(x = T.data$Sample,
                        levels = as.factor(str_split(SHcluster.Results$Sample[SHcluster.Results$order], pattern = ":", simplify = T)[, 2]),
                        ordered = TRUE)
#T.sample <- str_split(unique(T.data$Sample), pattern = ":", simplify = T)[, 2]
#T.sample <- paste0(formatC(c(1:6, length(T.sample):7), width = 2, flag = "0"), T.sample)
#T.data$Sample <- rep(T.sample, each = length(unique(T.data$MCLNum)))

T.data$FDR[T.data$FDR > 5e-2] <- 400
T.data$FDR[T.data$FDR <= 5e-2 & T.data$FDR > 5e-5] <- 300
T.data$FDR[T.data$FDR <= 5e-5 & T.data$FDR > 5e-10] <- 200
T.data$FDR[T.data$FDR <= 5e-10 & T.data$FDR > 5e-20] <- 100
T.data$FDR[T.data$FDR <= 5e-20] <- 0
#Create heatmap plot
heatmap.plot <- ggplot(data = T.data, aes(x = MCLNum, y = Sample, fill = FDR))
heatmap.plot <- heatmap.plot +geom_tile(color = "black", size = 0.1)
heatmap.plot <- heatmap.plot + scale_fill_gradient(low="seagreen",high="white")
heatmap.plot <- heatmap.plot + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
heatmap.plot <- heatmap.plot + theme(axis.ticks.y = element_blank())
heatmap.plot <- heatmap.plot + theme(legend.position = "top")
heatmap.plot <- heatmap.plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
heatmap.plot <- heatmap.plot + theme(legend.title = element_blank())
heatmap.plot <- heatmap.plot + theme(legend.text = element_text(size=7))
heatmap.plot <- heatmap.plot + theme(legend.position = 'none')
#check:VennDiagram--------------
T.data <- CY15_EnrichemntOfGenesList %>% select(-starts_with("pvalue"), -contains("BTH"), -contains("MCLNum"), -contains("CY16"), -contains("CY20"))
publish <- rownames(T.data)[which(apply(T.data[, 1:5] < 0.05, MARGIN = 1, FUN = sum) != 0)]
CY15 <- rownames(T.data)[which(apply(T.data[, 7:ncol(T.data)] < 0.05, MARGIN = 1, FUN = sum) != 0)]
T.Venn <- list(publish, CY15)
names(T.Venn) <- c("published", "CY15")
CY.Venn <- venn(T.Venn)
names(attr(CY.Venn,"intersections")) <- paste0(names(attr(CY.Venn,"intersections")), ",")
T.data <- data.frame(CY = str_split(names(unlist(attr(CY.Venn,"intersections"))), pattern = ",", simplify = T)[, 1],
                     MCLNum = unlist(attr(CY.Venn,"intersections"), use.names = F),
                     stringsAsFactors = F)
#save----
ggsave("~/bigdata/yasue/PCCOfPCC/CY/Image/EnrichmentOfDEGs/CY15andGenesList_heatmap.png", heatmap.plot, width = 24, height = 16)
write.table(T.data, file = "~/bigdata/yasue/PCCOfPCC/CY/Table/PublishedAndCY15_VennDiagram.txt", sep = "\t", quote = F, row.names = F)
ggsave(filename = "~/bigdata/yasue/PCCOfPCC/CY/Image/EnrichmentOfDEGs/PublishedAndCY15_VennDiagram.png", plot = plot(CY.Venn), width = 6, height = 4)
#remove object----
rm(list = ls())