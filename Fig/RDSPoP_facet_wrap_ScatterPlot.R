"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/Fig/RDSPoP_facet_wrap_ScatterPlot.R"
library(ggplot2)
library(stringr)
path <- list.files("~/bigdata/yasue/PCCOfPCC/RDS/PoP_facet_wrap/")
path <- paste0("~/bigdata/yasue/PCCOfPCC/RDS/PoP_facet_wrap/", path)
T.MCLNum <- str_split(path, pattern = "CY15MCLNum", simplify = T)[, 2]
T.MCLNum <- str_split(T.MCLNum, pattern = "CY15", simplify = T)[, 1]
j <- 1
for(j in j:length(T.MCLNum)){
  PoP <- readRDS(path[j])
  g <- ggplot(PoP, aes(x=CY15, y=OpenData))
  g <- g + geom_point()
  g <- g + xlab("CY15 PCC") + ylab("OpenData PCC")
  g <- g + facet_wrap(~ sample, ncol = 4)
  g <- g + xlim(-1, 1) + ylim(-1, 1)
  g <- g + theme_classic()
  #g <- g + stat_density2d(aes(alpha=..density..), geom = "tile",contour = FALSE)
  g <- g + stat_density2d(aes(fill = ..level..), alpha = 0.3, geom = "polygon")
  g <- g + scale_fill_continuous(low = "grey", high = "red", space = "Lab", name = "g = 0")
  g <- g +  scale_colour_discrete(guide = FALSE)
  g <- g + theme(strip.text = element_text(size = 9))
  g <- g + theme(axis.text.x = element_text(size = 10))
  g <- g + theme(axis.text.y = element_text(size = 10))
  plot(g)
  ImageTite <- paste0("~/bigdata/yasue/PCCOfPCC/CY/Image/PoP_ScatterPlot/facet_wrap/CY15MCLNum", T.MCLNum[j], "CY15PCCvsOpenDataPCC.png")
  ggsave(filename = ImageTite, plot = g, width = 12, height = 6)
}
