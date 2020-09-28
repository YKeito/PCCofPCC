#"/home/yasue/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/DAPSeq2/ATTEEDII_PeakIndexPhenotype.R"
#load package----
library(dplyr)
library(stringr)
#FioreDB Input, processing----
#Individual project
I.data <- read.table("/home/yasue/bigdata/yasue/FioreDB/IndividualProject.txt", sep = "\t", header = T, stringsAsFactors = F)
I.data$Num.of.lines[is.na(I.data$Num.of.lines)] <- 0
I.data$Num.of.lines.showing.the.phenotype[is.na(I.data$Num.of.lines.showing.the.phenotype)] <- 0
I.data$Sum.of.intesity.within.same.ID[is.na(I.data$Sum.of.intesity.within.same.ID)] <- 0
I.data$Gene <- toupper(I.data$Gene)
I.data <- I.data %>% filter(Num.of.lines.showing.the.phenotype > 1)
#Bulk project
B.data <- read.table("/home/yasue/bigdata/yasue/FioreDB/BulkProject.txt", sep = "\t", header = T, stringsAsFactors = F)
B.data[is.na(B.data$Phenotype.code), 4:7] <- B.data[is.na(B.data$Phenotype.code), 5:8]
B.data <- B.data %>% select(-X)
B.data$Causal.gene <- toupper(B.data$Causal.gene)
B.data <- B.data[grep("AT", B.data$Causal.gene), ]
PhenotypeData <- data.frame(AGI = c(I.data$Gene, B.data$Causal.gene),
                            PhynotypeCode = c(I.data$Phenotype.code, B.data$Phenotype.code),
                            Tissue = c(I.data$Tissue, B.data$Tissue),
                            Phynotype = c(I.data$Phenotype, B.data$Phenotype),
                            stringsAsFactors = F
)
#join PeakIndex, FioreDB----
T.data <- readRDS(file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/PeakIndex_thver.rds")
PeakIndex.Phenotype <- T.data %>% mutate(sAGI.Phenotype = rep(FALSE, nrow(T.data)), tAGI.Phenotype = rep(FALSE, nrow(T.data)))
T.AGI <- unique(PhenotypeData$AGI)
i <- 1
for(i in i:length(T.AGI)){
  PeakIndex.Phenotype$sAGI.Phenotype[T.AGI[i] == PeakIndex.Phenotype$sAGI] <- TRUE
  PeakIndex.Phenotype$tAGI.Phenotype[T.AGI[i] == PeakIndex.Phenotype$tAGI] <- TRUE
  i <- i+1
}
saveRDS(object = PeakIndex.Phenotype, file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/20190809DAPSeq_GenomeWide/PeakIndexPhenotype.rds")
write.table(x = PeakIndex.Phenotype, file = "~/PeakIndexPhenotype.txt", sep = "\t", quote = F, row.names = F)
PeakIndex.Phenotype %>% subset(PeakIndex.Phenotype %>% select(contains("Phenotype")) %>% apply(MARGIN = 1, FUN = sum) != 0) %>% select(sAGI, tAGI) %>% unlist(use.names = FALSE) %>% unique() %>% length()