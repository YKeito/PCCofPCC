#"/home/yasue/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/FioreDB/PeakIndexPhenotype_FioreDB_MatchRate.R"
#input data----
TAIR10 <- read.table("/home/yasue/bigdata/yasue/TAIR10/TAIR10_functional_descriptions", sep = "\t", header = TRUE, quote = "", fill = TRUE, stringsAsFactors = FALSE)
PeakIndex.Phenotype <- readRDS("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/PeakIndexPhenotype.rds")
sort.matrix <- readRDS(file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/FioreDB/UMAPCluster_FioreDB_SortMatrix.rds")
#formatting data----
TAIR10$Model_name <- paste0(str_sub(TAIR10$Model_name, start = 1, end = 9), "_", str_sub(TAIR10$Model_name, start = 11))
TAIR10 <- TAIR10 %>% filter(grepl("_1", Model_name))
TAIR10$Model_name <- str_sub(TAIR10$Model_name, start = 1, end = 9)
T.data <- PeakIndex.Phenotype %>% filter(sAGI.Phenotype == TRUE & tAGI.Phenotype == TRUE)
#processing data----
Match.Rate <- c()
T.Phenotype <- c()
T.Tissue <- c()
T.ID <- c()
i <- 1
for(i in i:nrow(T.data)){
  predata <- sort.matrix[T.data$sAGI[i] == sort.matrix$AGI | T.data$tAGI[i] == sort.matrix$AGI, 3:ncol(sort.matrix)]
  T.Match <- predata %>% lapply(FUN = as.numeric) %>% data.frame() %>% apply(MARGIN = 2, FUN = sum) == 2
  Match.Rate <- c(Match.Rate, sum(T.Match)/ncol(predata))
  T.ID <- c(T.ID, paste("NA", colnames(sort.matrix)[c(2+which(T.Match))], collapse = " | "))
  T.Phenotype <- c(T.Phenotype, paste("NA", sort.matrix[1, c(2+which(T.Match))], collapse = " | "))
  T.Tissue <- c(T.Tissue, paste("NA", sort.matrix[2, c(2+which(T.Match))], collapse = " | "))
}
T.ID[nchar(T.ID) != 3] <- gsub("NA ", "", T.ID[nchar(T.ID) != 3])
T.Phenotype[nchar(T.Phenotype) != 3] <- gsub("NA ", "", T.Phenotype[nchar(T.Phenotype) != 3])
T.Tissue[nchar(T.Tissue) != 3] <- gsub("NA ", "", T.Tissue[nchar(T.Tissue) != 3])
T.data <- T.data %>% mutate(Match.Rate = Match.Rate,
                            FioreDBID = T.ID,
                            Tissue = T.Tissue,
                            Phenotype = T.Phenotype)
saveRDS(object = T.data, file = "~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/FioreDB/PeakIndexPhenotype_FioreDB_MatchRate.rds")
write.table(T.data, file = "~/bigdata/yasue/PCCOfPCC/ATTEDII/Table/FioreDB/ColocalizedPhenotype/PeakIndexPhenotype_FioreDB_MatchRate.txt", sep = "\t", quote = F, row.names = F)