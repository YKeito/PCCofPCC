#"~/Nakano_RNAseq/network_analysis/script/PCCofPCC/Fig/DeterminationOfThreshold.R"
before <- proc.time()
#processing data----------
T.CY <- c("CY15", "CY16", "CY20")
n <- 1
for(n in n:length(T.CY)){
  title <- paste0("~/bigdata/yasue/PCCOfPCC/RDS/", T.CY[n], "_PCC_dataframe.rds")
  CY15.PCC <- readRDS(file = title)
  ####number of edge####
  #base edges 3,565,785
  Numedges <- c()
  value <- c(5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5)
  total_i <- length(value)
  #th:FDR & possitive correlation
  i <- 1
  for(i in i:total_i){
    Numedges <- c(Numedges, nrow(CY15.PCC[CY15.PCC$q_value < value[i] & CY15.PCC$interaction_value > 0, ]))
    i <- i+1
  }
  png(paste0("~/bigdata/yasue/PCCOfPCC/CY/Image/DeterminatinOfThresholds/", T.CY[n], "_FDR_Possitive.png"), width = 549, height = 437)
  plot(value, Numedges,
       main=paste0(T.CY[n], " : Number of edges in various thresholds of PCC"),
       xlab="FDR & possitive correlation used for cut-off", ylab="Numbers of edges"
  )
  abline(h=6.0e+05, lty=2)
  dev.off()
  #th:FDR
  Numedges <- c()
  value <- c(5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5)
  total_i <- length(value)
  i <- 1
  for(i in i:total_i){
    Numedges <- c(Numedges, nrow(CY15.PCC[CY15.PCC$q_value < value[i], ]))
    i <- i+1
  }
  png(paste0("~/bigdata/yasue/PCCOfPCC/CY/Image/DeterminatinOfThresholds/", T.CY[n], "_FDR.png"), width = 549, height = 437)
  plot(value, Numedges,
       main=paste0(T.CY[n], " : Number of edges in various thresholds of PCC"),
       xlab="FDR used for cut-off", ylab="Numbers of edges"
  )
  abline(h=6.0e+05, lty=2)
  dev.off()
  
  #th:possitive correlation
  Numedges <- c()
  value <- c(0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99)
  total_i <- length(value)
  i <- 1
  for(i in i:total_i){
    Numedges <- c(Numedges, nrow(CY15.PCC[CY15.PCC$interaction_value >= value[i], ]))
    i <- i+1
  }
  png(paste0("~/bigdata/yasue/PCCOfPCC/CY/Image/DeterminatinOfThresholds/", T.CY[n], "_possitive_PCC.png"), width = 549, height = 437)
  plot(value, Numedges,
       main=paste0(T.CY[n], " : Number of edges in various thresholds of PCC"),
       xlab="PCC used for cut-off", ylab="Numbers of edges"
  )
  abline(h=6.0e+05, lty=2)
  dev.off()
  
  #th:abs PCC
  Numedges <- c()
  value <- c(0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99)
  total_i <- length(value)
  i <- 1
  for(i in i:total_i){
    Numedges <- c(Numedges, nrow(CY15.PCC[abs(CY15.PCC$interaction_value) >= value[i], ]))
    i <- i+1
  }
  png(paste0("~/bigdata/yasue/PCCOfPCC/CY/Image/DeterminatinOfThresholds/", T.CY[n], "_absPCC.png"), width = 549, height = 437)
  plot(value, Numedges,
       main=paste0(T.CY[n], " : Number of edges in various thresholds of PCC"),
       xlab="|PCC| used for cut-off", ylab="Numbers of edges"
  )
  abline(h=6.0e+05, lty=2)
  dev.off()
  print(n)
}
#elapsed time------------------------------------------
after <- proc.time()
print(after - before)#2465.704 sec,41.09507 min
#remove object-------
rm(list = ls())