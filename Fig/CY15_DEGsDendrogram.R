"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/Fig/DEGs_Dendrogram.R"
#package----
library(dplyr)
library(tidyr)
library(factoextra)
library(ggplot2)
library(stringr)
#input data----
CY15_EnrichemntOfGenesList <- readRDS("~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_EnrichemntOfGenesList.rds")
#processing data----
T.data <- CY15_EnrichemntOfGenesList %>% select(-starts_with("pvalue"), -contains("BTH"), -contains("MCLNum"), -contains("CY16"), -contains("CY20"))
T.data <- T.data[apply(T.data < 0.05, MARGIN = 1, FUN = sum) != 0, ]
#Hclust----
#MCLNum
res.hc <- eclust(x = T.data,
                 "hclust",
                 k = ncol(T.data),
                 method = "euclidean",
                 graph = FALSE
)
#Sample
Sample.hc <- eclust(x = t(T.data),
                 "hclust",
                 k = 2,
                 method = "euclidean",
                 graph = FALSE
)
#Dendrogram----
#MCLNum
g <- fviz_dend(res.hc,
               cex = 0.5,
               color_labels_by_k = TRUE,
               show_labels = TRUE,
               ggtheme = theme_classic(),
               horiz = FALSE,
               #rect = TRUE,
               rect_fill = TRUE,
               type = "rectangle",
               main = NULL
)
g <- g + theme(axis.text=element_text(size=20))
g <- g + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
T.HCL <- data.frame(MCLNum = names(res.hc$cluster),
                    Hcluster = res.hc$cluster,
                    order = res.hc$order,
                    row.names = NULL)
#Sample
g.sample <- fviz_dend(Sample.hc,
                      cex = 0.5,
                      color_labels_by_k = TRUE,
                      show_labels = TRUE,
                      ggtheme = theme_classic(),
                      horiz = FALSE,
                      #rect = TRUE,
                      rect_fill = TRUE,
                      type = "rectangle",
                      main = NULL
)
g.sample <- g.sample + theme(axis.text=element_text(size=10))
g.sample <- g.sample + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
Sample.HCL <- data.frame(Sample = names(Sample.hc$cluster),
                         Hcluster = Sample.hc$cluster,
                         order = Sample.hc$order,
                         row.names = NULL)
#option----
#k:クラスター数の指定
#k_colors, palette:色指定
#show_labels:横軸の名前を表示するか
#color_labels_by_k:クラスターのグループに応じて色を付けるか
#horiz縦軸と横軸を反転するかしないか
#rect:クラスターを囲む
#rect_fill:rectで囲んだ範囲を色塗する
#type:図の形式の指定
#lwd:dendrogramの線の太さ変更
#label_cols:下の文字の色指定
#save----
#MCLNum
ggsave("~/bigdata/yasue/PCCOfPCC/CY/Image/EnrichmentOfDEGs/CY15_MCLNumDendrogram.png", g, width = 20, height = 10)
saveRDS(T.HCL, "~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_MCLNumDendrogram.rds")
write.table(T.HCL, "~/bigdata/yasue/PCCOfPCC/CY/Table/CY15_MCLNumDendrogram.txt", sep = "\t", append = F, quote = F, row.names = F, col.names = T)
#Sample
ggsave("~/bigdata/yasue/PCCOfPCC/CY/Image/EnrichmentOfDEGs/CY15_SampleDendrogram.png", g.sample, width = 20, height = 10)
saveRDS(Sample.HCL, "~/bigdata/yasue/PCCOfPCC/RDS/Table/CY15_SampleDendrogram.rds")
write.table(Sample.HCL, "~/bigdata/yasue/PCCOfPCC/CY/Table/CY15_SampleDendrogram.txt", sep = "\t", append = F, quote = F, row.names = F, col.names = T)