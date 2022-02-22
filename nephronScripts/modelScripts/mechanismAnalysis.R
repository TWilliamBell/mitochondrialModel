if (isFALSE(grepl("modelScripts", getwd()))) {
  setwd("nephronScripts/modelScripts")
}

files <- c(paste0('../results/resultsmTALMech', 1:6, '.csv'), 
           paste0('../results/resultsmTALMech', 1:4, 'C3.csv'))

mechanism <- lapply(files, data.table::fread)
mechanism <- lapply(mechanism, function(x) tail(x, 1))
mechanism <- sapply(mechanism, function(x) x$ATP_c)

baseline <- tail(read.csv("../results/resultsmTAL.csv"), 1)$ATP_c
comparison <- mechanism/baseline

hypoxia <- round(comparison[1:6], 3)
hypoxiamTAL <- min(read.csv("../results/resultsHypoxiaTAL0.csv")$ATP_c)
#hypoxiamTAL/baseline
comparisonPTHypox <- hypoxiaPT/baseline
## PT Hypoxia:
## 0.19
hypoxiaPT <- tail(read.csv("../results/resultsExtremeHypoxia.csv")$ATP_c, 1)

C3 <- round(comparison[7:10], 3)
C3PT <- read.csv("../results/tailsMitDisPT.csv")
C3PT <- C3PT[C3PT$nums, ]
iterProd <- read.csv("../results/iterProd.csv")
iterProd$nums <- 1:256
C3PT <- dplyr::inner_join(iterProd, C3PT, by = "nums")
C3PT <- C3PT[C3PT$`0` == 1 & C3PT$`1` == 0.25 & C3PT$`2` == 1 & C3PT$`3` == 1, ]$ATP_c
comparisonPTC3 <- C3PT/baseline
## PT C3:
## 0.69

