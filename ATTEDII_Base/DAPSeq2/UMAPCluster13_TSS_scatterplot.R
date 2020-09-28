library(ggplot2)
library(dplyr)
#"/home/yasue/Nakano_RNAseq/network_analysis/script/PCCofPCC/ATTEDII_Base/DAPSeq2/UMAPCluster13_TSS_scatterplot.R"
diff.TSS <- readRDS(file = "/home/yasue/bigdata/yasue/PCCOfPCC/ATTEDII/RDS/Table/UMAPCluster13/UMAPCluster13_741genes_TSS.rds")
T.sample <- diff.TSS$TFName %>% unique()
T.sample <- sort(T.sample)
h <- 1
for(h in h:2){
  T.value <- diff.TSS %>% filter(TFName == T.sample[h]) %>% select(value)
  i <- 3
  for(i in i:length(T.sample)){
    df <- cbind(x = T.value,
                y = diff.TSS %>% filter(TFName == T.sample[i]) %>% select(value)
    )
    colnames(df) <- c("A", "B")
    df <- df %>% filter(A <= 1000 & A >= -1000, B <= 1000 & B >= -1000)
    a <- lm(formula = A~B, data = df)
    a <- summary(a)
    g <- ggplot(df, aes(x = A, y = B))
    g <- g + geom_point()
    g <- g + theme_bw()
    g <- g + xlab(T.sample[h]) + ylab(T.sample[i])
    g <- g + ggtitle(paste0(T.sample[h], ":", T.sample[i]))
    g <- g + ggtitle(paste0(T.sample[h], ":", T.sample[i], ", R^2 = ", formatC(a$r.squared, digits = 2)))
    g <- g + stat_smooth(method = "lm", se = FALSE, colour = "red", size = 1)
    #plot(g)
    title <- paste0("~/bigdata/yasue/UMAP_Project/Image/TSS_scatterplot/UMAPCluster13/", T.sample[h], "_", T.sample[i], ".png")
    ggsave(filename = title, plot = g)
  }
  h <- h+1
}
