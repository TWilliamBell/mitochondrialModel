mitDis <- read.csv("../results/tailsBroadSim.csv")

normalTail <- tail(read.csv("../results/resultsATP.csv"), 1)

iterProd <- read.csv("../results/iterprod.csv")
joinedTable <- dplyr::inner_join(iterProd, mitDis, by = "nums")
joinedTable <- joinedTable[joinedTable$PO2 == 1.0 &
                             joinedTable$Glycolysis == 0.00034, ]
mitDis <- joinedTable

Q <- c(mitDis$QH2_x, normalTail$QH2_x)/6.49e-3
C <- c(mitDis$Cred_i, normalTail$Cred_i)/4.39e-3
N <- c(mitDis$NADH_x, normalTail$NADH_x)/9.4e-3

hist(N, main = "", xlab = "Reduced Proportion of NADH/NAD+ Pool",
     xlim = c(0, 1), col = "peachpuff1")

hist(Q, main = "", xlab = "Reduced Proportion of Coenzyme Q Pool",
     xlim = c(0, 1), col = "peachpuff1")
abline(v = normalTail$QH2_x/6.49e-3, col = "red")

#pdf("../dataVis/cytCOXPHOS.pdf")
hist(C, main = "", xlab = "Reduced Proportion of Cytochrome C Pool",
     xlim = c(0, 1), col = "peachpuff1")
abline(v = normalTail$Cred_i/4.39e-3, col = "red")
#dev.off()
