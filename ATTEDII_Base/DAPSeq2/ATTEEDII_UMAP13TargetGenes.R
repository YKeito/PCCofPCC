#"/home/yasue/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/DAPSeq2/ATTEEDII_UMAP13TargetGenes.R"
PeakIndex.Phenotype <- readRDS("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/PeakIndexPhenotype.rds")
PeakIndex.Phenotype <- PeakIndex.Phenotype %>% filter(Cluster == "UMAPCluster13")
DAPSeq.TFName <- readRDS("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/DAPSeq.TFName.rds")
T.family <- str_split(DAPSeq.TFName$Category, pattern = "/", simplify = T)[, 1]
names(T.family) <- str_split(DAPSeq.TFName$Category, pattern = "/", simplify = T)[, 2]
Colocalized.target <- c()
DAPSeq.target <- c()
g <- 1
for(g in g:length(PeakIndex.Phenotype)){
  title <-   paste0("/home/yasue/bigdata/yasue/DAPSeq/dap_data_v4/genes/", 
                    str_split(DAPSeq.TFName$Category[grep(PeakIndex.Phenotype$source[g], DAPSeq.TFName$Category)], pattern = "narrowPeak", simplify = T)[, 1],
                    "nS_targets.txt")
  source.data <- read.table(title, sep = "\t", header = T)
  title <-   paste0("/home/yasue/bigdata/yasue/DAPSeq/dap_data_v4/genes/", 
                    str_split(DAPSeq.TFName$Category[grep(PeakIndex.Phenotype$target[g], DAPSeq.TFName$Category)], pattern = "narrowPeak", simplify = T)[, 1],
                    "nS_targets.txt")
  target.data <- read.table(title, sep = "\t", header = T)
  T.target <- intersect(source.data$target.at_id, target.data$target.at_id)
  Colocalized.target <- rbind(Colocalized.target,
                              data.frame(sourceTF = rep(PeakIndex.Phenotype$source[g], times = length(T.target)),
                                         sourceTF.AGI = rep(PeakIndex.Phenotype$sAGI[g], times = length(T.target)),
                                         targetTF = rep(PeakIndex.Phenotype$target[g], times = length(T.target)),
                                         targetTF.AGI = rep(PeakIndex.Phenotype$tAGI[g], times = length(T.target)),
                                         Colocalized.target = T.target,
                                         stringsAsFactors = F)
  )
  h <- 1
  for(h in h:2){
    if(h == 1){
      DAPSeq.target <- rbind(DAPSeq.target,
                             data.frame(TF = PeakIndex.Phenotype$source[g],
                                        AGI = PeakIndex.Phenotype$sAGI[g],
                                        target.AGI = setdiff(source.data$target.at_id, target.data$target.at_id),
                                        stringsAsFactors = F
                             )
      )
    }else{
      DAPSeq.target <- rbind(DAPSeq.target,
                             data.frame(TF = PeakIndex.Phenotype$target[g],
                                        AGI = PeakIndex.Phenotype$tAGI[g],
                                        target.AGI = setdiff(target.data$target.at_id, source.data$target.at_id),
                                        stringsAsFactors = F
                             )
      )
    }
    h <- h+1
  }
  g <- g+1
}
DAPSeq.target <- DAPSeq.target[!duplicated(DAPSeq.target), ]
saveRDS(object = DAPSeq.target, file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/DAPSeq_targetgenes/UMAPCLuster13_NotOverlapTargetGenes.rds")
saveRDS(object = Colocalized.target, file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/DAPSeq_targetgenes/UMAPCLuster13_OverlapTargetGenes.rds")
write.table(x = Colocalized.target, file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/Table/DAPSeq/ColocalizedTarget.txt", sep = "\t")

