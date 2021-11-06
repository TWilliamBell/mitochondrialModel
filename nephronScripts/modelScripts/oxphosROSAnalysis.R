drugSimParam <- as.data.frame(feather::read_feather(
  "./results/iterProdDrugSim.feather"))

colnames(drugSimParam) <- c("k", "CI", "CIII", "CIV", "ATPSynthase",
                            "XHle", "PDHEq")

drugSimTails <- read.csv("./results/tailsDrugSim.csv")
normal <- data.table::fread("./results/resultsATP.csv")
normalTail <- tail(normal, 1)
normalTail$V1 <- normalTail$t <- NULL
normalTail <- unlist(normalTail)
drugSimTails$X <- drugSimTails$V1 <- drugSimTails$t <- NULL


colnames(drugSimTails)[1] <- "k"
joinedTable <- dplyr::inner_join(drugSimParam, drugSimTails, by = 'k')
drugSimParam <- drugSimParam[!(drugSimParam$PDHEq == 1
                               & drugSimParam$XHle == 1), ]

Q <- joinedTable$QH2_x/2.418e-3
C <- joinedTable$Cred_i/1.956e-3

hist(Q, main = "", xlab = "Reduced Proportion of Coenzyme Q",
     xlim = c(0, 1))
abline(v = normalTail[10]/2.418e-3, col = "red")

hist(C, main = "", xlab = "Reduced Proportion of Cytochrome C",
     xlim = c(0, 1))
abline(v = normalTail[27]/1.956e-3, col = "red")
