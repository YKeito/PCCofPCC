"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/Fig/DEGs_heatmap.R"
#package----
library(ggplot2)
library(stringr)
#input data----
CY15_EnrichemntOfGenesList <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_EnrichemntOfGenesList.rds")
#processing data----
T.data <- CY15_EnrichemntOfGenesList[CY15_EnrichemntOfGenesList$`qvalue:CY15_12h_FDR0.05` < 0.05 | CY15_EnrichemntOfGenesList$`qvalue:CY15_1h_FDR0.05` < 0.05 |
                                       CY15_EnrichemntOfGenesList$`qvalue:CY15_24h_FDR0.05` < 0.05 | CY15_EnrichemntOfGenesList$`qvalue:CY15_3h_FDR0.05` < 0.05, ]
T.data <- T.data[T.data$MCLNum < 115, ]
T.data <- T.data %>% select(starts_with("qvalue:CY15"))
T.data <- data.frame(MCLNum = rep(c(paste0("00", rownames(T.data)[nchar(rownames(T.data)) == 1]), 
                                    paste0("0", rownames(T.data)[nchar(rownames(T.data)) == 2]),
                                    rownames(T.data)[nchar(rownames(T.data)) == 3]),
                                  times = ncol(T.data)),
                     value = unlist(T.data),
                     sample = rep(paste0("0", c(ncol(T.data):1), str_split(colnames(T.data), pattern = ":", simplify = T)[, 2]), each = nrow(T.data)))
T.data$value[T.data$value > 5e-2] <- 500
T.data$value[T.data$value <= 5e-2 & T.data$value > 5e-5] <- 400
T.data$value[T.data$value <= 5e-5 & T.data$value > 5e-10] <- 300
T.data$value[T.data$value <= 5e-10 & T.data$value > 5e-20] <- 200
T.data$value[T.data$value <= 5e-20 & T.data$value > 5e-40] <- 100
T.data$value[T.data$value <= 5e-40] <- 0
#ggplot----
g <- ggplot(data = T.data, aes(x = MCLNum, y = sample, fill = value))
g <- g + geom_tile(color = "black", size = 0.5)
g <- g + scale_fill_gradient(low="seagreen",high="white")
g <- g + theme_linedraw()
g <- g + theme_set(theme_bw(base_size = 20))
g <- g +　theme(axis.text=element_text(size=24))
g <- g + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
g <- g + theme(legend.position="top")
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot(g)
#save data----
ggsave("~/bigdata/yasue/PCCOfPCC/CY/Image/EnrichmentOfDEGs/CY15DEGs_heatmap.png", g, width = 24, height = 8)
#remove object----
rm(list = ls())