#"/home/yasue/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/FioreDB/PhynotypeList.R"
#install package----
#install.packages("dplyr")
#install.packages("stringr")
#load package----
library(dplyr)
library(stringr)
#Individual project----
I.data <- read.table("/home/yasue/bigdata/yasue/FioreDB/IndividualProject.txt", sep = "\t", header = T, stringsAsFactors = F)
I.data$Num.of.lines[is.na(I.data$Num.of.lines)] <- 0
I.data$Num.of.lines.showing.the.phenotype[is.na(I.data$Num.of.lines.showing.the.phenotype)] <- 0
I.data$Sum.of.intesity.within.same.ID[is.na(I.data$Sum.of.intesity.within.same.ID)] <- 0
I.data$Gene <- toupper(I.data$Gene)
I.data <- I.data %>% filter(Num.of.lines.showing.the.phenotype > 1)
#Bulk project----
B.data <- read.table("/home/yasue/bigdata/yasue/FioreDB/BulkProject.txt", sep = "\t", header = T, stringsAsFactors = F)
B.data[is.na(B.data$Phenotype.code), 4:7] <- B.data[is.na(B.data$Phenotype.code), 5:8]
B.data <- B.data %>% select(-X)
B.data$Causal.gene <- toupper(B.data$Causal.gene)
B.data <- B.data[grep("AT", B.data$Causal.gene), ]
PhynotypeData <- data.frame(AGI = c(I.data$Gene, B.data$Causal.gene), 
                            PhynotypeCode = c(I.data$Phenotype.code, B.data$Phenotype.code),
                            Phynotype = c(I.data$Phenotype, B.data$Phenotype),
                            stringsAsFactors = F
                            )
#Phynotype List----
T.ID <- unique(PhynotypeData$PhynotypeCode)
PhynotypeList <- list()
h <- 1
for(h in h:length(T.ID)){
  T.data <- PhynotypeData %>% filter(PhynotypeCode == T.ID[h])
  if(nrow(T.data) != 0){
    PhynotypeList <- c(PhynotypeList, list(T.data$AGI %>% unique()))
    names(PhynotypeList)[length(PhynotypeList)] <- paste0(T.ID[h], ":", unique(T.data$Phynotype))
  }
  h <- h+1
}
#save data----
saveRDS(object = PhynotypeList, file = "/home/yasue/bigdata/yasue/FioreDB/FioreDB/PhynotypeList.rds")
