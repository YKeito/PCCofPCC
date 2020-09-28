"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/Fig/Supple_NumberOfNodesInMCLNum.R"
MasterTable <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_MCLNumNodeTable.rds")
T.MCLNum <- MasterTable %>% filter(MCLNum < 115)
T.MCLNum <- sort(T.MCLNum$MCLNum)

T.data <- data.frame(MCLNum = T.MCLNum,
                     stringsAsFactors = F)

g <- ggplot(T.data, aes(x = MCLNum))
g <- g + geom_histogram(binwidth = 1)
g <- g + theme(axis.text=element_text(size=14))
g <- g + theme(axis.title=element_text(size=14))
g <- g + theme(axis.title.x = element_blank())
g <- g + theme(axis.title.y = element_blank())

ggsave("~/bigdata/yasue/PCCOfPCC/CY/Image/NetworkFig/NumberOfNodesInMCLNum.png", plot = g, width = 5, height = 6)