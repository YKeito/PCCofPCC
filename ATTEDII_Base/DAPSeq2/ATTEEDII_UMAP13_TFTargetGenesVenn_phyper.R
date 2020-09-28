"/home/yasue/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/DAPSeq2/ATTEEDII_UMAP13_TFTargetGenesVenn_phyper.R"
#CY.data----
CY.data <- read.table("~/Nakano_RNAseq/network_analysis/base/CY_base/allRNASeq_union.txt", sep = "\t", header = T)
CY.data <- CY.data %>% select(ends_with("q_value"), -contains("48h"))
#target----
T.matrix <- readRDS(file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/UMAPCluster13_DAPSeqTargetTable.rds")
#data processing----
HY5.target <- T.matrix %>% filter(AT5G11260 == 1) %>% select("AGI")
HY5.target <- HY5.target$AGI
AT3G42860.target <- T.matrix %>% filter(AT3G42860 == 1) %>% select("AGI")
AT3G42860.target <- AT3G42860.target$AGI
WRKYs <- T.matrix[, match(setdiff(colnames(T.matrix)[nchar(colnames(T.matrix)) == 9], c("AT5G11260", "AT3G42860", "AT5G09410")), colnames(T.matrix))] == 1
WRKYs <- WRKYs %>% apply(MARGIN = 1, FUN = sum)
WRKYs.target <- T.matrix$AGI[which(WRKYs != 0)]
data <- list(HY5 = HY5.target, zfGRF = AT3G42860.target, WRKYs = WRKYs.target)
Venn <- venn(data)
names(attr(Venn,"intersections")) <- paste0(names(attr(Venn,"intersections")), ",")
T.Venn <- data.frame(AGI = unlist(attr(Venn, "intersections")),
                     Sample = str_split(names(unlist(attr(Venn, "intersections"))), pattern = ",", simplify = T)[, 1], 
                     stringsAsFactors = F,
                     row.names = NULL)
T.sample <- T.Venn$Sample %>% unique()
N <- nrow(T.matrix)
CY.EnrichmentOfTFTargets <- c()
i <- 1
for(i in i:length(T.sample)){
  target.data <- T.Venn %>% filter(Sample == T.sample[i])
  T.AGI <- target.data$AGI
  n <- T.AGI %>% length()
  j <- 1
  for(j in j:ncol(CY.data)){
    CY.DEGs <- CY.data %>% select(colnames(CY.data)[j])< 0.05
    CY.DEGs <- rownames(CY.data)[CY.DEGs]
    M <- CY.DEGs %>% length()
    x <- intersect(T.AGI, CY.DEGs) %>% length()
    CY.EnrichmentOfTFTargets <- rbind(CY.EnrichmentOfTFTargets, 
                                      data.frame(TF = T.sample[i],
                                                 CY = colnames(CY.data)[j],
                                                 pvalue = phyper(x-1, M, N-M, n, lower.tail = F),
                                                 x = x,
                                                 M = M,
                                                 N = N,
                                                 n = n,
                                                 stringsAsFactors = F
                                                 )
                                      )
    j <- j+1
    if(j == c(ncol(CY.data)+1)){
      CY.DEGs <- rownames(CY.data)
      M <- CY.DEGs %>% length()
      x <- intersect(T.AGI, CY.DEGs) %>% length()
      CY.EnrichmentOfTFTargets <- rbind(CY.EnrichmentOfTFTargets, 
                                        data.frame(TF = T.sample[i],
                                                   CY = "allCY",
                                                   pvalue = phyper(x-1, M, N-M, n, lower.tail = F),
                                                   x = x,
                                                   M = M,
                                                   N = N,
                                                   n = n,
                                                   stringsAsFactors = F
                                        )
      )
    }
  }
}
saveRDS(object = T.Venn, file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/UMAPCluster13/UMAPCluster13_DAPSeqTargetVennTable.rds")
saveRDS(object = CY.EnrichmentOfTFTargets, file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/UMAPCluster13/UMAPCluster13_DAPSeqTargetVennEnrichmentTable.rds")