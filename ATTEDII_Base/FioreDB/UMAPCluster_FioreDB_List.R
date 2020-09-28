#"/home/yasue/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/FioreDB/UMAPCluster_FioreDB_List.R"
TAIR10 <- read.table("/home/yasue/bigdata/yasue/TAIR10/TAIR10_functional_descriptions", sep = "\t", header = TRUE, quote = "", fill = TRUE, stringsAsFactors = FALSE)
TAIR10$Model_name <- paste0(str_sub(TAIR10$Model_name, start = 1, end = 9), "_", str_sub(TAIR10$Model_name, start = 11))
TAIR10 <- TAIR10 %>% filter(grepl("_1", Model_name))
TAIR10$Model_name <- str_sub(TAIR10$Model_name, start = 1, end = 9)
UMAPCluster <- "UMAPCluster5"
PeakIndex.Phenotype <- readRDS("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/PeakIndexPhenotype.rds")
PeakIndex.Phenotype <- PeakIndex.Phenotype %>% filter(sAGI.Phenotype == TRUE & tAGI.Phenotype == TRUE, Cluster == UMAPCluster)
sort.matrix <- readRDS(file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/FioreDB/UMAPCluster_FioreDB_SortMatrix.rds")
sort.matrix <- sort.matrix %>% filter(UMAPCluster == UMAPCluster | UMAPCluster == "")
T.AGI <- union(PeakIndex.Phenotype$sAGI, PeakIndex.Phenotype$tAGI)
Colocalized.Phenotype.data <- c()
g <- 1
for(g in g:length(T.AGI)){
  T.PeakIndex.Phenotype <- PeakIndex.Phenotype %>% filter(tAGI == T.AGI[g] | sAGI == T.AGI[g])
  T.Colocalized.Phenotype.data <- c()
  h <- 1
  for(h in h:nrow(T.PeakIndex.Phenotype)){
    T.data <- sort.matrix %>% filter(AGI == T.PeakIndex.Phenotype$sAGI[h] | AGI == T.PeakIndex.Phenotype$tAGI[h])
    if(nrow(T.data) != 0){
      T.data <- T.data[, 3:ncol(T.data)]
      count <- T.data %>% lapply(FUN = as.numeric) %>% data.frame() %>% apply(MARGIN = 2, FUN = sum) == 2
      if(sum(count) != 0){
        T.Colocalized.Phenotype.data <- rbind(T.Colocalized.Phenotype.data, 
                                              sort.matrix %>% filter(AGI == T.PeakIndex.Phenotype$sAGI[h] | AGI == T.PeakIndex.Phenotype$tAGI[h])
                                              )
      }
    }
    h <- h+1
  }
  if(!is.null(T.Colocalized.Phenotype.data)){
    T.Colocalized.Phenotype.data <- T.Colocalized.Phenotype.data %>% distinct(AGI, .keep_all = TRUE)
    T.logic <- T.Colocalized.Phenotype.data[, 3:ncol(T.Colocalized.Phenotype.data)] %>% lapply(FUN = as.numeric) %>% data.frame() %>% apply(MARGIN = 2, FUN = sum) >= 2
    T.Colocalized.Phenotype.data <- rbind(sort.matrix[1:2, ], T.Colocalized.Phenotype.data)
    T.Colocalized.Phenotype.data$UMAPCluster[1] <- T.AGI[g]
    T.Colocalized.Phenotype.data <- T.Colocalized.Phenotype.data[, c(1:2, which(T.logic)+2)]
    T.annotation <- TAIR10$Short_description[match(T.Colocalized.Phenotype.data$AGI, TAIR10$Model_name)]
    T.annotation[c(1, 2)] <- ""
    T.Colocalized.Phenotype.data <- T.Colocalized.Phenotype.data %>% mutate(annotation = T.annotation)
    T.Colocalized.Phenotype.data <- T.Colocalized.Phenotype.data[, c(1, ncol(T.Colocalized.Phenotype.data), 2:c(ncol(T.Colocalized.Phenotype.data)-1))]

    title <- paste0("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/Table/FioreDB/ColocalizedPhenotype/", UMAPCluster, "/UMAPCluster_FioreDB_Colocalized_", T.AGI[g], ".txt")
    write.table(T.Colocalized.Phenotype.data, file = title, sep = "\t", quote = F, row.names = F)
  }
  g <- g+1
}
