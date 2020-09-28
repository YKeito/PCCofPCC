#"/home/yasue/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/UMAP/UMAPCluster_STRING_List.R"
TAIR10 <- read.table("/home/yasue/bigdata/yasue/TAIR10/TAIR10_functional_descriptions", sep = "\t", header = TRUE, quote = "", fill = TRUE, stringsAsFactors = FALSE)
TAIR10$Model_name <- paste0(str_sub(TAIR10$Model_name, start = 1, end = 9), "_", str_sub(TAIR10$Model_name, start = 11))
TAIR10 <- TAIR10 %>% filter(grepl("_1", Model_name))
TAIR10$Model_name <- str_sub(TAIR10$Model_name, start = 1, end = 9)
STRING <- read.table("/home/yasue/bigdata/yasue/STRING/3702.protein.actions.v11.0.txt", sep = "\t", header = TRUE, quote = "", fill = TRUE, stringsAsFactors = FALSE)
STRING$item_id_a <- paste0(str_sub(STRING$item_id_a, start = 6, end = 14), "_", str_sub(STRING$item_id_a, start = 16))
STRING$item_id_b <- paste0(str_sub(STRING$item_id_b, start = 6, end = 14), "_", str_sub(STRING$item_id_b, start = 16))
STRING <- STRING %>% filter(grepl("_1", item_id_a))
STRING <- STRING %>% filter(grepl("_1", item_id_b))
STRING$item_id_a <- str_sub(STRING$item_id_a, start = 1, end = 9)
STRING$item_id_b <- str_sub(STRING$item_id_b, start = 1, end = 9)
DAPSeq.TFName <- readRDS("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/DAPSeq.TFName.rds")
temp <- intersect(match(DAPSeq.TFName$AGI, STRING$item_id_a), match(DAPSeq.TFName$AGI, STRING$item_id_b))
temp <- temp[!is.na(temp)]
T.STRING <- STRING[temp, ]
T.STRING <- T.STRING  %>% mutate(UMAP_a = DAPSeq.TFName$UMAPCluster[match(T.STRING$item_id_a, DAPSeq.TFName$AGI)],
                                 UMAP_b = DAPSeq.TFName$UMAPCluster[match(T.STRING$item_id_b, DAPSeq.TFName$AGI)],
                                 annotation_a = TAIR10$Short_description[match(T.STRING$item_id_a, TAIR10$Model_name)],
                                 annotation_b = TAIR10$Short_description[match(T.STRING$item_id_b, TAIR10$Model_name)])
write.table(T.STRING, file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/Table/STRING/UMAPCluster_STRING.txt", sep = "\t", quote = F, row.names = F)
