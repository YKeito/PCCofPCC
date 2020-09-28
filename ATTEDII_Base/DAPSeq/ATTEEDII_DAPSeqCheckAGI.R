"/home/yasue/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/DAPSeq/ATTEEDII_DAPSeqCheckAGI.R"
DAPSeq.TFAGI <- readRDS("/home/yasue/bigdata/yasue/DAPSeq/DAPSeq_TF.AGI.rds") #DAP-Seqの解析対象TFは387遺伝子
T.data <- DAPSeq.TFAGI %>% filter(AGI == "AT1G27730" | AGI == "AT1G28370" | AGI == "AT2G23320" | AGI == "AT3G15210" | AGI == "AT4G01250" | AGI == "AT4G18880" | AGI == "AT4G31800")
temp <- combn(x = T.data$pwd, m = 2)
h <- 1 #137, UMAPCluster4
for(h in h:ncol(temp)){
  base.pwd <- paste0("/home/yasue/bigdata/yasue/DAPSeq/dap_data_v4/peaks/", temp[1, h])
  base.peak <- read.table(base.pwd)
  colnames(base.peak) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
  target.pwd <- paste0("/home/yasue/bigdata/yasue/DAPSeq/dap_data_v4/peaks/", temp[2, h])
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
      print(paste0("UMAPCluster", g, ", sample:", ncol(temp)-h, " chrom:", i, "NumRows:", length(T.locus)-j))
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
  title <- paste0("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/DAPSeq_genomewide/CheckList/DAPSeq_base", str_split(temp[1, h], pattern = "/", simplify = T)[, 2], "&target", str_split(temp[2, h], pattern = "/", simplify = T)[, 2])
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
  title <- paste0("/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/Image/DAPSeq_genomewide/CheckList/DAPSeq_base", str_split(temp[1, h], pattern = "/", simplify = T)[, 2], "&target", str_split(temp[2, h], pattern = "/", simplify = T)[, 2])
  ggsave(paste0(title, ".png"), df, width = 6, height = 4)
  h <- h+1
}