#setwd("../results")

mitDis <- data.table::fread("../results/tailsMitDismTALSim.csv")
mitDis$V1 <- mitDis$t <- NULL
colnames(mitDis)[1] <- "k"

iterProd <- feather::read_feather('../results/iterProd.feather')
iterProd <- cbind(1:256, iterProd, rep(1, 256), rep(0, 256), 
                  rep(1, 256))
colnames(iterProd) <- c("k", "CI", "CIII", "CIV", "ATPSynthase",
                        "XHle", "glycolysisLevel", "O2Level")
mitDis <- dplyr::inner_join(iterProd, mitDis)

## Identify univariate cases
univar <- rep(FALSE, nrow(mitDis))
j <- 0
for (i in mitDis$k) {
  j <- j+1
  if (prod(iterProd[i, 2:5]) %in% c(0.25, 0.5, 0.75, 1) & sum(iterProd[i, 2:5]) > 3) {
    univar[j] <- TRUE
  }
}

## Consider ATP concentration
ATP <- mitDis$ATP_c*1000
Complexes <- as.matrix(mitDis[ , 2:5])

pdf("../dataVis/univarPlatesmTALATP.pdf")
par(mfrow = c(2, 2), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)

plot(c(0.25, 0.5, 0.75, 1), ATP[(mitDis$CI != 1.0 | 
                                   matrixStats::rowProds(Complexes) == 1.)
                                & univar],
     ylim = c(0, max(ATP, na.rm = T)), ylab = "Cytosolic ATP (mM)",
     xlim = c(0.2, 1.), xlab = "Fraction of Typical Activity",
     col = "blue", main = "Complex I")
CI <- splinefun(c(0.25, 0.5, 0.75, 1), 
                ATP[(mitDis$CI != 1.0 | 
                       matrixStats::rowProds(Complexes) == 1.)
                    & univar])
curve(CI, from = 0.2, to = 1, add = T, col = "blue")
abline(a = min(ATP), b = 0, col = "red")

plot(c(0.25, 0.5, 0.75, 1), ATP[(mitDis$CIII != 1.0 | 
                                   matrixStats::rowProds(Complexes) == 1.)
                                & univar],
     ylim = c(0, max(ATP, na.rm = T)), ylab = "Cytosolic ATP (mM)",
     xlim = c(0.2, 1.), xlab = "Fraction of Typical Activity", 
     col = "blue", main = "Complex III")
#rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4)-0.0001, c(0.25, 0.5, 0.75, 1)+0.05, 
#     ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 2] != 1.0) & univar],
#     col = "springgreen")
CIII <- splinefun(c(0.25, 0.5, 0.75, 1), ATP[(mitDis$CIII != 1.0 | 
                                                matrixStats::rowProds(Complexes) == 1.)
                                             & univar])
curve(CIII, from = 0.2, to = 1, add = T, col = "blue")
abline(a = min(ATP), b = 0, col = "red")

plot(c(0.25, 0.5, 0.75, 1), ATP[(mitDis$CIV != 1.0 | 
                                   matrixStats::rowProds(Complexes) == 1.)
                                & univar],
     ylim = c(0, max(ATP, na.rm = T)), ylab = "Cytosolic ATP (mM)",
     xlim = c(0.2, 1.), xlab = "Fraction of Typical Activity", 
     col = "blue", main = "Complex IV")
CIV <- splinefun(c(0.25, 0.5, 0.75, 1), ATP[(mitDis$CIV != 1.0 | 
                                               matrixStats::rowProds(Complexes) == 1.)
                                            & univar])
curve(CIV, from = 0.2, to = 1., add = T, col = "blue")
abline(a = min(ATP), b = 0, col = "red")

plot(c(0.25, 0.5, 0.75, 1), ATP[(mitDis$CI != 1.0 | 
                                   matrixStats::rowProds(Complexes) == 1.)
                                & univar],
     ylim = c(0, max(ATP, na.rm = T)), ylab = "Cytosolic ATP (mM)",
     xlim = c(0.2, 1.), xlab = "Fraction of Typical Activity",
     col = "blue", main = "ATP Synthase")
#rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4)-0.0001, c(0.25, 0.5, 0.75, 1)+0.05, 
#     ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 4] != 1.0) & univar],
#     col = "springgreen")
ATPcurve <- splinefun(c(0.25, 0.5, 0.75, 1), ATP[(mitDis$CI != 1.0 | 
                                                    matrixStats::rowProds(Complexes) == 1.)
                                                 & univar])
curve(ATPcurve, from = 0.2, to = 1, add = T, col = "blue")
abline(a = ATP[which.min(ATP)], b = 0, col = "red")

dev.off()
