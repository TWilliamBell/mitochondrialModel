if (isFALSE(grepl("modelscripts", getwd()))) {
  setwd("modelscripts")
}

files <- c(paste0('../results/resultsmTALMech', 1:6, '.feather'), 
           paste0('../results/resultsmTALMech', 1:4, 'C3.feather'))

mechanism <- lapply(files, feather::read_feather)
mechanism <- lapply(mechanism, function(x) tail(x, 1))
mechanism <- sapply(mechanism, function(x) x$ATP_c)

baseline <- tail(feather::read_feather("../results/resultsmTAL.feather"), 1)$ATP_c
comparison <- mechanism/baseline

hypoxia <- round(comparison[1:6], 3)
hypoxiaPT <- min(read.csv("../results/resultsHypoxia0.csv")$ATP_c)
hypoxiamTAL <- min(read.csv("../results/resultsHypoxiaTAL0.csv")$ATP_c)
#hypoxiamTAL/baseline
comparisonPTHypox <- hypoxiaPT/baseline
## PT Hypoxia:
## 0.19
## Hypoxia:
## 0.73 k_O2, 0.36 J_AtC, 0.98 V_mito, 0.22 V_mito & J_AtC,
## 0 J_AtC & k_O2, 0.48 V_mito & k_O2

C3 <- round(comparison[7:10], 3)
C3PT <- read.csv("../results/tailsMitDisPT.csv")
C3PT <- C3PT[C3PT$nums, ]
iterProd <- as.data.frame(feather::read_feather("../results/iterProd.feather"))
iterProd$nums <- 1:256
C3PT <- dplyr::inner_join(iterProd, C3PT, by = "nums")
C3PT <- C3PT[C3PT$`0` == 1 & C3PT$`1` == 0.25 & C3PT$`2` == 1 & C3PT$`3` == 1, ]$ATP_c
comparisonPTC3 <- C3PT/baseline
## PT C3:
## 0.69
## C3:
## 0.99 k_O2, 0.91 J_AtC, 0.96 V_mito, 0.60 J_AtC & V_mito, 0.53 all three
