if (!grepl("nephronScripts/modelScripts", getwd())) {
  setwd("nephronScripts/modelScripts")
}

## In order to run this file, you should first run simulateMitDis1.py and then collectMitDisPT.R &
## writeProduct.py

mitDis <- list()

trueVec <- rep(TRUE, 256)

for (i in 1:256) {
  if (paste0("resultsDiseaseATP", i, ".csv") %in% dir("../results")) {
    mitDis[[i]] <- data.table::fread(paste0("../results/resultsDiseaseATP", i, ".csv"))
  } else {
    trueVec[i] <- FALSE
  }
}

par(cex.lab = 1.5)

iterProd <- as.data.frame(read.csv("../results/iterProd.csv"))
iterProd$X <- NULL

trueInd <- (1:256)[trueVec]

tails <- list()
conc <- names(unlist(tail(mitDis[[5]], 1))[-(1:2)]) ## Chose the first one where results
## were found

for (i in trueInd) {
  tails[[i]] <- tail(mitDis[[i]], 1)#unname(unlist(tail(mitDis[[i]], 1))[-(1:2)])
}

for (i in setdiff(1:256, trueInd)) {
  tails[[i]] <- data.frame(matrix(nrow = 1, ncol = 65))
  colnames(tails[[i]]) <- colnames(tails[[250]])
}
## Identify univariate cases

univar <- rep(FALSE, 256)

for (i in trueInd) {
  if (prod(iterProd[i, ]) %in% c(0.25, 0.5, 0.75, 1) & sum(iterProd[i, ]) > 3) {
    univar[i] <- TRUE
  }
}

## Consider ATP concentration
ATP <- sapply(tails, function(x) if (is.null(x)) {return(NA)} else {unlist(x)[37]})*1000

# ## Complex I
# pdf("../dataVis/complexIATPunivar.pdf")
# par(cex.lab = 1.5)
# plot(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
#                                   iterProd[ , 1] != 1.0) & univar],
#      ylim = c(0, max(ATP, na.rm = T)), ylab = "ATP Concentration",
#      xlim = c(0.2, 1.1), xlab = "Fraction of Typical Activity",
#      cex = 0.1, main = "Complex I")
# rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
#      ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |iterProd[ , 1] != 1.0) & univar],
#      col = "springgreen")
# abline(a = ATP[which.min(ATP)], b = 0, col = "red")
# dev.off()
# 
# ## Complex III
# pdf("../dataVis/complexIIIATPunivar.pdf")
# par(cex.lab = 1.5)
# plot(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
#                                    iterProd[ , 2] != 1.0) & univar],
#      ylim = c(0, max(ATP, na.rm = T)), ylab = "ATP Concentration",
#      xlim = c(0.1, 1.1), xlab = "Fraction of Typical Activity", 
#      cex = 0.1,
#      main = "Complex III")
# rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
#      ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 2] != 1.0) & univar],
#      col = "springgreen")
# abline(a = ATP[which.min(ATP)], b = 0, col = "red")
# dev.off()
# 
# ## Complex IV
# pdf("../dataVis/complexIVATPunivar.pdf")
# par(cex.lab = 1.5)
# plot(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
#                                    iterProd[ , 3] != 1.0) & univar],
#      ylim = c(0, max(ATP, na.rm = T)), ylab = "ATP Concentration",
#      xlim = c(0.1, 1.1), xlab = "Fraction of Typical Activity", cex = 0.1,
#      main = "Complex IV")
# rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
#      ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 3] != 1.0) & univar],
#      col = "springgreen")
# abline(a = ATP[which.min(ATP)], b = 0, col = "red")
# dev.off()
# 
# ## ATP Synthase
# pdf("../dataVis/atpsynthaseATPunivar.pdf")
# par(cex.lab = 1.5)
# plot(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
#                                    iterProd[ , 4] != 1.0) & univar],
#      ylim = c(0, max(ATP, na.rm = T)), ylab = "ATP Concentration",
#      xlim = c(0.1, 1.1), xlab = "Fraction of Typical Activity",
#      cex = 0.1, main = "ATP Synthase")
# rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
#      ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 4] != 1.0) & univar],
#      col = "springgreen")
# abline(a = ATP[which.min(ATP)], b = 0, col = "red")
# dev.off()

## Figure 3.9

pdf("../dataVis/univarPlatePTATP.pdf")
par(cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)

plot(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                   iterProd[ , 1] != 1.0) & univar], col = "green", ylim = c(0, 2.5),
     xlim = c(0.2, 1), xlab = "Fraction of Typical Activity", ylab = "Cytosolic ATP (mM)")
CI <- splinefun(c(0.25, 0.5, 0.75, 1), 
                ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                       iterProd[ , 1] != 1.0) & univar])
curve(CI, from = 0.2, to = 1, add = T, col = "green")

points(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                   iterProd[ , 2] != 1.0) & univar], col = 'blue')
#rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4)-0.0001, c(0.25, 0.5, 0.75, 1)+0.05, 
#     ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 2] != 1.0) & univar],
#     col = "springgreen")
CIII <- splinefun(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                                iterProd[ , 2] != 1.0) & univar])
curve(CIII, from = 0.2, to = 1, add = T, col = "blue")

points(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                   iterProd[ , 3] != 1.0) & univar], col = "purple")
#rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4)-0.0001, c(0.25, 0.5, 0.75, 1)+0.05, 
#     ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 3] != 1.0) & univar],
#     col = "springgreen")
CIV <- splinefun(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                               iterProd[ , 3] != 1.0) & univar])
curve(CIV, from = 0.2, to = 1., add = T, col = "purple")

points(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                   iterProd[ , 4] != 1.0) & univar], col = "orange")
#rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4)-0.0001, c(0.25, 0.5, 0.75, 1)+0.05, 
#     ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 4] != 1.0) & univar],
#     col = "springgreen")
ATPcurve <- splinefun(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                                    iterProd[ , 4] != 1.0) & univar])
curve(ATPcurve, from = 0.2, to = 1, add = T, col = "orange")

abline(h = min(ATP), col = "red")

legend(x = 0.65, y = 0.75, fill = c("orange", "green", "blue", "purple"),
       legend = c("ATP Synthase", "Complex I", "Complex III", "Complex IV"), cex = 1.5)
dev.off()

pdf("../dataVis/univarPlates.pdf")
par(mfrow = c(2, 2), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)

plot(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                   iterProd[ , 1] != 1.0) & univar],
     ylim = c(0, max(ATP, na.rm = T)), ylab = "Cytosolic ATP (mM)",
     xlim = c(0.2, 1.), xlab = "Fraction of Typical Activity",
     col = "blue", main = "Complex I")
CI <- splinefun(c(0.25, 0.5, 0.75, 1), 
                ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                              iterProd[ , 1] != 1.0) & univar])
curve(CI, from = 0.2, to = 1, add = T, col = "blue")
#rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4)-0.0001, c(0.25, 0.5, 0.75, 1)+0.05, 
#     ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |iterProd[ , 1] != 1.0) & univar],
#     col = "springgreen")
abline(a = ATP[which.min(ATP)], b = 0, col = "red")

plot(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                   iterProd[ , 2] != 1.0) & univar],
     ylim = c(0, max(ATP, na.rm = T)), ylab = "Cytosolic ATP (mM)",
     xlim = c(0.2, 1.), xlab = "Fraction of Typical Activity", 
     col = "blue", main = "Complex III")
#rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4)-0.0001, c(0.25, 0.5, 0.75, 1)+0.05, 
#     ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 2] != 1.0) & univar],
#     col = "springgreen")
CIII <- splinefun(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                                iterProd[ , 2] != 1.0) & univar])
curve(CIII, from = 0.2, to = 1, add = T, col = "blue")
abline(a = ATP[which.min(ATP)], b = 0, col = "red")

plot(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                   iterProd[ , 3] != 1.0) & univar],
     ylim = c(0, max(ATP, na.rm = T)), ylab = "Cytosolic ATP (mM)",
     xlim = c(0.2, 1.), xlab = "Fraction of Typical Activity", 
     col = "blue", main = "Complex IV")
#rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4)-0.0001, c(0.25, 0.5, 0.75, 1)+0.05, 
#     ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 3] != 1.0) & univar],
#     col = "springgreen")
CIV <- splinefun(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                               iterProd[ , 3] != 1.0) & univar])
curve(CIV, from = 0.2, to = 1., add = T, col = "blue")
abline(a = ATP[which.min(ATP)], b = 0, col = "red")

plot(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                   iterProd[ , 4] != 1.0) & univar],
     ylim = c(0, max(ATP, na.rm = T)), ylab = "Cytosolic ATP (mM)",
     xlim = c(0.2, 1.), xlab = "Fraction of Typical Activity",
     col = "blue", main = "ATP Synthase")
#rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4)-0.0001, c(0.25, 0.5, 0.75, 1)+0.05, 
#     ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 4] != 1.0) & univar],
#     col = "springgreen")
ATPcurve <- splinefun(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                                    iterProd[ , 4] != 1.0) & univar])
curve(ATPcurve, from = 0.2, to = 1, add = T, col = "blue")
abline(a = ATP[which.min(ATP)], b = 0, col = "red")

dev.off()

par(mfrow = c(1, 1))

## Consider the potential gradient
dPsi <- sapply(tails, function(x) if (is.null(x)) {return(NA)} else {x[[4]]})

## Complex I
pdf("../dataVis/complexIdPsiunivar.pdf")
par(cex.lab = 1.5)
plot(c(0.25, 0.5, 0.75, 1), dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                   iterProd[ , 1] != 1.0) & univar],
     ylim = c(150, max(dPsi, na.rm = T)), ylab = "Electrical Potential Gradient",
     xlim = c(0.1, 1.1), xlab = "Fraction of Typical Activity",
     cex = 0.1, main = "Complex I")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(149.57, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
     dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 1] != 1.0) & univar],
     col = "springgreen")
abline(a = dPsi[which.min(dPsi)], b = 0, col = "red")
dev.off()

## Ifosfamide

## Figure 3.11
pdf("../dataVis/ifosfamide.pdf")
par(cex.lab = 1.5)
plot(c(0.5, 1), dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
            iterProd[ , 1] == 0.5) & univar],
     ylim = c(155, max(dPsi, na.rm = T)), ylab = "Electrical Potential Gradient",
     xlim = c(0.1, 1.1), xlab = "Fraction of Typical Activity",
     cex = 0.1, main = "Ifosfamide")
rect(c(0.5, 1)-0.05, rep(154.77, 2), c(0.5, 1)+0.05, 
     dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 1] == 0.5) & univar],
     col = "springgreen")
dev.off()

## Complex III
pdf("../dataVis/complexIIIdPsiunivar.pdf")
par(cex.lab = 1.5)
plot(c(0.25, 0.5, 0.75, 1), dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                   iterProd[ , 2] != 1.0) & univar],
     ylim = c(150, max(dPsi, na.rm = T)), ylab = "Electrical Potential Gradient",
     xlim = c(0.1, 1.1), xlab = "Fraction of Typical Activity", cex = 0.1,
     main = "Complex III")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(149.57, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
     dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 2] != 1.0) & univar],
     col = "springgreen")
abline(a = dPsi[which.min(dPsi)], b = 0, col = "red")
dev.off()
## Complex IV

pdf("../dataVis/complexIVdPsiunivar.pdf")
par(cex.lab = 1.5)
plot(c(0.25, 0.5, 0.75, 1), dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                   iterProd[ , 3] != 1.0) & univar],
     ylim = c(150, max(dPsi, na.rm = T)), ylab = "Electrical Potential Gradient",
     xlim = c(0.1, 1.1), xlab = "Fraction of Typical Activity", cex = 0.1,
     main = "Complex IV")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(149.57, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
     dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 3] != 1.0) & univar],
     col = "springgreen")
abline(a = dPsi[which.min(dPsi)], b = 0, col = "red")
dev.off()
## ATP Synthase

pdf("../dataVis/atpsynthasedPsiunivar.pdf")
par(cex.lab = 1.5)
plot(c(0.25, 0.5, 0.75, 1), dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                   iterProd[ , 4] != 1.0) & univar],
     ylim = c(150, max(dPsi, na.rm = T)), ylab = "Electrical Potential Gradient",
     xlim = c(0.1, 1.1), xlab = "Fraction of Typical Activity",
     cex = 0.1, main = "ATP Synthase")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(149.57, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
     dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 4] != 1.0) & univar],
     col = "springgreen")
abline(a = dPsi[which.min(dPsi)], b = 0, col = "red")
dev.off()

pdf("../dataVis/univarPlatesdPsi.pdf")
par(mfrow = c(2,2), cex.lab = 1.5)

plot(c(0.25, 0.5, 0.75, 1), dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                    iterProd[ , 1] != 1.0) & univar],
     ylim = c(150, max(dPsi, na.rm = T)), ylab = "Electrical Potential Gradient",
     xlim = c(0.1, 1.1), xlab = "Fraction of Typical Activity",
     cex = 0.1, col = "blue", main = "Complex I")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(149.57, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
     dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 1] != 1.0) & univar],
     col = "springgreen")
abline(a = dPsi[which.min(dPsi)], b = 0, col = "red")

plot(c(0.25, 0.5, 0.75, 1), dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                    iterProd[ , 2] != 1.0) & univar],
     ylim = c(150, max(dPsi, na.rm = T)), ylab = "Electrical Potential Gradient",
     xlim = c(0.1, 1.1), xlab = "Fraction of Typical Activity", cex = 0.1,
     col = "blue", main = "Complex III")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(149.57, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
     dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 2] != 1.0) & univar],
     col = "springgreen")
abline(a = dPsi[which.min(dPsi)], b = 0, col = "red")

plot(c(0.25, 0.5, 0.75, 1), dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                    iterProd[ , 3] != 1.0) & univar],
     ylim = c(150, max(dPsi, na.rm = T)), ylab = "Electrical Potential Gradient",
     xlim = c(0.1, 1.1), xlab = "Fraction of Typical Activity", cex = 0.1,
     col = "blue", main = "Complex IV")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(149.57, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
     dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 3] != 1.0) & univar],
     col = "springgreen")
abline(a = dPsi[which.min(dPsi)], b = 0, col = "red")

plot(c(0.25, 0.5, 0.75, 1), dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                    iterProd[ , 4] != 1.0) & univar],
     ylim = c(150, max(dPsi, na.rm = T)), ylab = "Electrical Potential Gradient",
     xlim = c(0.1, 1.1), xlab = "Fraction of Typical Activity",
     cex = 0.1, col = "blue", main = "ATP Synthase")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(149.57, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
     dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 4] != 1.0) & univar],
     col = "springgreen")
abline(a = dPsi[which.min(dPsi)], b = 0, col = "red")

dev.off()

## Look at scale of interactions
allResults <- do.call(rbind, tails)
allResults <- cbind(iterProd, allResults)
colnames(allResults)[1:4] <- c("CI", "CIII", "CIV", "ATPSynthase")

## The ATP response is well-approximated by a linear surface in CI, CIII, CIV, ATP Synthase, adding mixed terms doesn't
## improve the fit, suggesting that multiple causes interact additively

## Table 3.7 uses these values
LMATP <- lm(allResults$ATP_c ~ allResults$CI + allResults$CIII + allResults$CIV + allResults$ATPSynthase, allResults)
print(summary(LMATP))

LMATPMix <- lm(allResults$ATP_c ~ allResults$CI + allResults$CIII + allResults$CIV + allResults$CI*allResults$CIII*allResults$CIV*allResults$ATPSynthase, allResults)
print(summary(LMATPMix))

library(ggplot2)
container <- allResults
container$CIII <- as.factor(container$CIII)

## Figure 3.10
pdf("../dataVis/atpComplexIIImultivar.pdf")
p <- ggplot(container, aes(x = ATP_c*1000)) +
  geom_histogram(aes(fill = CIII), bins = 24) +
  xlab("Cytosolic ATP (mM)") +
  ylab("Frequency") +
  theme(legend.text = element_text(size = 18), axis.title = element_text(size = 18),
        legend.title = element_text(size = 18), axis.text = element_text(size = 18)) +
  xlim(c(0, 2.5))
p$labels$fill <- "Relative \nComplex III \nActivity"
print(p)
dev.off()

container <- allResults[allResults$CIII == 1, ]
container$CIV <- as.factor(container$CIV)

pdf("../dataVis/atpComplexIVmultivar.pdf")
ggplot(container, aes(x = ATP_c*1000)) +
  geom_histogram(aes(color = CIV, fill = CIV), bins = 24) +
  xlab("ATP Concentration (mmol)") +
  ylab("Frequency") +
  theme(legend.text = element_text(size = 18), axis.title = element_text(size = 18),
        legend.title = element_text(size = 18), axis.text = element_text(size = 18))
dev.off()

## Table 3.8 uses these values
LMPSI <- lm(allResults$ATP_c ~ allResults$CI + allResults$CIII + allResults$CIV, allResults)
print(summary(LMPSI))

LMPSIMix <- lm(allResults$ATP_c ~ allResults$CI + allResults$CIII + allResults$CIV +
                 allResults$CI*allResults$CIII*allResults$CIV, allResults)
print(summary(LMPSIMix))

container <- allResults
container$CIII <- as.factor(container$CIII)

pdf("../dataVis/dPsiComplexIIImultivar.pdf")
ggplot(container, aes(x = dPsi)) +
  geom_histogram(aes(color = CIII, fill = CIII), bins = 25) +
  xlab("dPsi") +
  ylab("Frequency") +
  theme(legend.text = element_text(size = 18), axis.title = element_text(size = 18),
        legend.title = element_text(size = 18), axis.text = element_text(size = 18))
dev.off()

container <- allResults[allResults$CIII == 1, ]
container$CIV <- as.factor(container$CIV)

pdf("../dataVis/dPsiComplexIVmultivar.pdf")
ggplot(container, aes(x = dPsi)) +
  geom_histogram(aes(color = CIV, fill = CIV)) +
  xlab("dPsi") +
  ylab("Frequency") +
  theme(legend.text = element_text(size = 18), axis.title = element_text(size = 18),
        legend.title = element_text(size = 18), axis.text = element_text(size = 18))
dev.off()

