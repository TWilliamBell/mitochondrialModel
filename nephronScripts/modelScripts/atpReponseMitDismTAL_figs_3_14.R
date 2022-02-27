## In order to run this file, you should first run simulateMitDismTAL.py, and then collectMitDismTAL.R 
## & writeProduct.py

if (!grepl("nephronScripts/modelScripts", getwd())) {
  setwd("nephronScripts/modelScripts")
}

mitDis <- data.table::fread("../results/tailsMitDismTALSim.csv")
mitDis$V1 <- mitDis$t <- NULL
colnames(mitDis)[1] <- "k"

iterProd <- read.csv('../results/iterProd.csv')
iterProd$X <- NULL
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

## Figure 3.14

pdf("../dataVis/univarPlatemTALATP.pdf")
par(cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)

plot(c(0.25, 0.5, 0.75, 1), ATP[(mitDis$CI != 1.0 | 
                                   matrixStats::rowProds(Complexes) == 1.)
                                & univar],
     ylim = c(0, max(ATP, na.rm = T)), ylab = "Cytosolic ATP (mM)",
     xlim = c(0.2, 1.), xlab = "Fraction of Typical Activity",
     col = "green")
CI <- splinefun(c(0.25, 0.5, 0.75, 1), 
                ATP[(mitDis$CI != 1.0 | 
                       matrixStats::rowProds(Complexes) == 1.)
                    & univar])
curve(CI, from = 0.2, to = 1, add = T, col = "green")
abline(a = min(ATP), b = 0, col = "red")

points(c(0.25, 0.5, 0.75, 1), ATP[(mitDis$CIII != 1.0 | 
                                   matrixStats::rowProds(Complexes) == 1.)
                                & univar], col = "blue")
#rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4)-0.0001, c(0.25, 0.5, 0.75, 1)+0.05, 
#     ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 2] != 1.0) & univar],
#     col = "springgreen")
CIII <- splinefun(c(0.25, 0.5, 0.75, 1), ATP[(mitDis$CIII != 1.0 | 
                                                matrixStats::rowProds(Complexes) == 1.)
                                             & univar])
curve(CIII, from = 0.2, to = 1, add = T, col = "blue")

points(c(0.25, 0.5, 0.75, 1), ATP[(mitDis$CIV != 1.0 | 
                                   matrixStats::rowProds(Complexes) == 1.)
                                & univar], 
     col = "purple")
CIV <- splinefun(c(0.25, 0.5, 0.75, 1), ATP[(mitDis$CIV != 1.0 | 
                                               matrixStats::rowProds(Complexes) == 1.)
                                            & univar])
curve(CIV, from = 0.2, to = 1., add = T, col = "purple")

points(c(0.25, 0.5, 0.75, 1), ATP[(mitDis$ATPSynthase != 1.0 | 
                                   matrixStats::rowProds(Complexes) == 1.)
                                & univar], col = "orange")
#rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4)-0.0001, c(0.25, 0.5, 0.75, 1)+0.05, 
#     ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 4] != 1.0) & univar],
#     col = "springgreen")
ATPcurve <- splinefun(c(0.25, 0.5, 0.75, 1), ATP[(mitDis$ATPSynthase != 1.0 | 
                                                    matrixStats::rowProds(Complexes) == 1.)
                                                 & univar])
curve(ATPcurve, from = 0.2, to = 1, add = T, col = "orange")
legend(x = 0.65, y = 0.75, fill = c("orange", "green", "blue", "purple"),
       legend = c("ATP Synthase", "Complex I", "Complex III", "Complex IV"), cex = 1.5)
dev.off()


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

## For table 3.8

LMATP <- lm(mitDis$ATP_c ~ mitDis$CI + mitDis$CIII + mitDis$CIV, mitDis)
print(summary(LMATP))

LMATPMix <- lm(mitDis$ATP_c ~ mitDis$CI + mitDis$CIII + mitDis$CIV +
                 mitDis$CI*mitDis$CIII*mitDis$CIV, mitDis)
print(summary(LMATPMix))

LMPSI <- lm(mitDis$dPsi ~ mitDis$CI + mitDis$CIII + mitDis$CIV, mitDis)
print(summary(LMPSI))

LMPSIMix <- lm(mitDis$dPsi ~ mitDis$CI + mitDis$CIII + mitDis$CIV +
                 mitDis$CI*mitDis$CIII*mitDis$CIV, mitDis)
print(summary(LMPSIMix))

