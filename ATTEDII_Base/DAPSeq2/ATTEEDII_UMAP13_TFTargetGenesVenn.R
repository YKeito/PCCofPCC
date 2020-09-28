"/home/yasue/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/DAPSeq2/ATTEEDII_UMAP13_TFTargetGenesVenn.R"
#pacckage----
library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)
#processing data----
before <- proc.time()
#input data----
TAIR10 <- readRDS(file = "/home/yasue/bigdata/yasue/TAIR10_ShortName.rds")
DAPSeq.TFName <- readRDS("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/DAPSeq.TFName.rds")
PeakIndex <- readRDS(file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/PeakIndex.rds")
PeakIndex <- PeakIndex %>% select(-contains("PeakCount"))
T.data <- PeakIndex %>% 
  filter(Cluster == "UMAPCluster13") %>% 
  filter(PeakIndex > 2 & diffTFfamily == FALSE, PeakMax > 50)
T.data <- T.data %>% mutate(annotation.1 = rep(0, times = nrow(T.data)), annotation.2 = rep(0, times = nrow(T.data)))
i <- 1
for(i in i:nrow(T.data)){
  T.data$annotation.1[i] <- TAIR10$annotation[T.data$sAGI[i] == TAIR10$AGI]
  T.data$annotation.2[i] <- TAIR10$annotation[T.data$tAGI[i] == TAIR10$AGI]
  i <- i+1
}
T.matrix <- TAIR10 %>% cbind(matrix(rep(0, times = nrow(TAIR10)*length(union(T.data$sAGI, T.data$tAGI))),
                                    ncol = length(union(T.data$sAGI, T.data$tAGI))), deparse.level = )
colnames(T.matrix)[3:ncol(T.matrix)] <- union(T.data$sAGI, T.data$tAGI)
i <- 3
for(i in i:c(ncol(T.matrix)-1)){
  title <-   paste0("/home/yasue/bigdata/yasue/DAPSeq/dap_data_v4/genes/", 
                    str_split(DAPSeq.TFName$Category[match(colnames(T.matrix)[i], DAPSeq.TFName$AGI)], 
                              pattern = "narrowPeak", simplify = T)[, 1],
                    "nS_targets.txt")
  source.data <- read.table(title, sep = "\t", header = T, stringsAsFactors = F)
  T.match <- match(source.data$target.at_id, T.matrix$AGI)
  T.match <- T.match[!is.na(T.match)]
  T.matrix[T.match, i] <- 1
  i <- i+1
}
CY.data <- read.table("~/Nakano_RNAseq/network_analysis/base/CY_base/allRNASeq_union.txt", sep = "\t", header = T)
T.data <- CY.data %>% select(ends_with("q_value"), -contains("48h"))
T.data <- CY.data[which(apply(CY.data < 0.05, MARGIN = 1, FUN = sum) != 0), ]
CY.data <- T.data %>% select(1:36, -contains("48h"))
T.qvalue <- CY.data %>% select(ends_with("q_value"))
T.qvalue <- T.qvalue < 0.05
T.DEGs <- T.qvalue %>% apply(MARGIN = 2, FUN = which)
T.matrix <- T.matrix %>% cbind(matrix(rep(0, times = nrow(TAIR10)*length(T.DEGs)), 
                                      ncol = length(T.DEGs)))
colnames(T.matrix)[which(colnames(T.matrix) == "1"):ncol(T.matrix)] <- names(T.DEGs)
i <- 1
for(i in i:length(T.DEGs)){
  CY.DEGs <- rownames(CY.data)[T.DEGs[[i]]]
  T.matrix[match(CY.DEGs, T.matrix$AGI), names(T.DEGs)[i]] <- 1
  i <- i+1
}
saveRDS(object = T.matrix, file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/UMAPCluster13_DAPSeqTargetTable.rds")
#UMAPCluster13以外のWRKYの標的遺伝子群の重複----
HY5.target <- T.matrix %>% filter(AT5G11260 == 1) %>% select("AGI")
HY5.target <- HY5.target$AGI
AT3G42860.target <- T.matrix %>% filter(AT3G42860 == 1) %>% select("AGI")
AT3G42860.target <- AT3G42860.target$AGI
WRKYs <- T.matrix[, match(setdiff(colnames(T.matrix)[nchar(colnames(T.matrix)) == 9], c("AT5G11260", "AT3G42860", "AT5G09410")), colnames(T.matrix))] == 1
WRKYs <- WRKYs %>% apply(MARGIN = 1, FUN = sum)
WRKYs.target <- T.matrix$AGI[which(WRKYs != 0)]
data <- list(HY5 = HY5.target, zfGRF = AT3G42860.target, WRKYs = WRKYs.target)
Venn <- venn(data)
T.Venn <- data.frame(unlist(attr(Venn, "intersections")))

library(venneuler)
n.A <- length(HY5.target)
n.B <- length(AT3G42860.target)
n.C <- length(WRKYs.target)
n.AB <- intersect(HY5.target, AT3G42860.target) %>% length()
n.BC <- intersect(AT3G42860.target, WRKYs.target) %>% length()
n.CA <- intersect(WRKYs.target, HY5.target) %>% length()
n.ABC <- n.A + n.B + n.C - (n.AB + n.BC + n.CA)
v <- venneuler(c(HY5 = n.A, AT3G42860 = n.B, WRKYs = n.C, "HY5&AT3G42860" = n.AB, "AT3G42860&WRKYs" = n.BC, "WRKYs&HY5" = n.CA, "HY5&AT3G42860&WRKYs" = n.ABC))
plot(v, main = "Overlaps of ABC", col = c("darkgreen", "blue", "orange"))
#UMAPCluster13のWRKYの標的遺伝子群の重複----
WRKYs.target <- T.matrix[, match(setdiff(colnames(T.matrix)[nchar(colnames(T.matrix)) == 9], c("AT5G11260", "AT3G42860", "AT5G09410")), colnames(T.matrix))]
i <- 1
for(i in i:ncol(WRKYs.target)){
  temp <- T.matrix$AGI[WRKYs.target[, i] == 1]
  data <- list(HY5 = HY5.target, zfGRF = AT3G42860.target, WRKYs = temp)
  names(data)[3] <- TAIR10$annotation[match(colnames(WRKYs.target)[i], TAIR10$AGI)]
  title <- paste0("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/Image/UMAPCluster13_Fig/VennDiagram/HY5_AT3G42860_", names(data)[3], ".png")
  png(filename = title)
  venn(data)
  dev.off()
}
#UMAPCluster13以外のWRKYの標的遺伝子群の重複----
WRKY6.target <- read.table("/home/yasue/bigdata/yasue/DAPSeq/dap_data_v4/genes/WRKY_tnt/WRKY22_col/chr1-5/chr1-5_GEM_events.nS_targets.txt", sep = "\t", header = T)
WRKY6.target <- WRKY6.target$target.at_id
data <- list(HY5 = HY5.target, zfGRF = AT3G42860.target, WRKY6 = WRKY6.target)
Venn <- venn(data)
