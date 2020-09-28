"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/ATTEEDII_DAPSeqUMAP_Cluster.R"
before <- proc.time()
#input data----
UMAP.Cluster <- readRDS("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/20190809UMAPClusterList.rds")
T.data <- readRDS("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/EnrichmentAnalysis/ATTEDII_EnrichemntOfTFsTargetList.rds")
UMAP.TF <- intersect(UMAP.Cluster$Sample, T.data$TF) #UMAPの解析対象TFは570遺伝子
UMAP.TFAGI <- T.data$AGI[match(UMAP.TF, T.data$TF)]
names(UMAP.TFAGI) <- UMAP.TF
#DAP-Seq TFList----
DAPSeq.TFAGI <- readRDS("/home/yasue/bigdata/yasue/DAPSeq/DAPSeq_TF.AGI.rds") #DAP-Seqの解析対象TFは387遺伝子
T.AGI <- UMAP.TFAGI[match(DAPSeq.TFAGI$AGI, UMAP.TFAGI)]
T.AGI <- T.AGI[!is.na(T.AGI)] #UMAPの解析対象かつDAP-Seqの解析対象TFは316遺伝子
names(T.AGI) <- UMAP.Cluster$cluster[match(names(T.AGI), UMAP.Cluster$Sample)] #316遺伝子がUMAPのどのクラスターに分類されているか
DAPSeq.TFName <- data.frame(AGI = T.AGI,
                            DAPName = DAPSeq.TFAGI$TFName[match(T.AGI, DAPSeq.TFAGI$AGI)],
                            UMAPCluster = names(T.AGI),
                            Category = DAPSeq.TFAGI$pwd[match(T.AGI, DAPSeq.TFAGI$AGI)],
                            stringsAsFactors = F, row.names = NULL)
rm(list = ls()[ls() != "DAPSeq.TFName"&ls() != "before"])
saveRDS(DAPSeq.TFName, "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/DAPSeq.TFName.rds")
#processing data----
g <- 2
for(g in g:length(unique(DAPSeq.TFName$UMAPCluster))){
  T.data <- DAPSeq.TFName %>% filter(UMAPCluster == g)
  temp <- combn(x = T.data$Category, m = 2)
  h <- 1 #137, UMAPCluster4
  for(h in h:ncol(temp)){
    base.pwd <- paste0("~/bigdata/yasue/DAPSeq/dap_data_v4/peaks/", temp[1, h])
    base.peak <- read.table(base.pwd)
    colnames(base.peak) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
    target.pwd <- paste0("~/bigdata/yasue/DAPSeq/dap_data_v4/peaks/", temp[2, h])
    target.peak <- read.table(target.pwd)
    colnames(target.peak) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
    allchrom.LengthDifference <- c()
    chr <- unique(base.peak$chrom) %>% sort()
    i <- 1
    for(i in i:length(chr)){
      base.data <- base.peak %>% filter(chrom == chr[i])
      target.data <- target.peak %>% filter(chrom == chr[i])
      S.locus <- str_split(base.data$name, pattern = ":", simplify = T)[, 2] %>% as.numeric()
      T.locus <- str_split(target.data$name, pattern = ":", simplify = T)[, 2] %>% as.numeric()
      T.LengthDifference <- c()
      j <- 1
      for(j in j:length(T.locus)){
        T.LengthDifference <- c(T.LengthDifference, sweep(x = S.locus %>% as.matrix(), MARGIN = 1, STATS = T.locus[j], FUN = "-") %>% as.vector())
        print(paste0("UMAPCluster", g, ", sample:", ncol(temp)-h, " chrom:", i, " NumRows:", length(T.locus)-j))
        j <- j+1
      }
      allchrom.LengthDifference <- rbind(allchrom.LengthDifference, data.frame(chrom = rep(chr[i], times = length(T.LengthDifference)),
                                                                               value = T.LengthDifference,
                                                                               SRegion = rep(S.locus, times = length(T.locus)),
                                                                               TRegion = rep(T.locus, each = length(S.locus)),
                                                                               stringsAsFactors = F)
                                         )
      i <- i+1
    }
    #histogram----
    allchrom.LengthDifference <- allchrom.LengthDifference %>% filter(value <= 3000 & value >= -3000)
    title <- paste0("~/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/UMAPCluster", g, "/DAPSeq_base", str_split(temp[1, h], pattern = "/", simplify = T)[, 2], "&target", str_split(temp[2, h], pattern = "/", simplify = T)[, 2])
    saveRDS(allchrom.LengthDifference, file = paste0(title, ".rds"))
    #df <- ggplot(T.data, aes(x = value, fill = chrom))
    df <- ggplot(allchrom.LengthDifference, aes(x = value))
    df <- df + geom_histogram(position="identity", alpha = 0.6, binwidth = 50)
    df <- df + theme_bw()
    #df <- df + facet_grid(chrom ~ .)
    df <- df + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
    df <- df + theme(legend.position = 'none')
    df <- df + theme(axis.text=element_text(size=18))
    #df <- df + theme(strip.text.y = element_text(size=18))
    df <- df + scale_x_continuous(breaks = seq(-3000, 3000, 1000))
    title <- paste0("~/bigdata/yasue/PCCOfPCC/ATTEDII/Image/20190809DAPSeq_GenomeWide/UMAPCluster", g, "/DAPSeq_base", str_split(temp[1, h], pattern = "/", simplify = T)[, 2], "&target", str_split(temp[2, h], pattern = "/", simplify = T)[, 2])
    ggsave(paste0(title, ".png"), df, width = 6, height = 4)
    h <- h+1
  }
  g <- g+1
}
#elapsed time----
after <- proc.time()
print(after - before)#67.649 sec
rm(list = ls())