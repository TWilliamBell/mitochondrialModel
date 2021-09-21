#setwd("./results")

mitDis <- list()

trueVec <- rep(TRUE, 256)

for (i in 1:256) {
  if (paste0("resultsDiseaseATP", i, ".csv") %in% dir()) {
    mitDis[[i]] <- data.table::fread(paste0("resultsDiseaseATP", i, ".csv"))
  } else {
    trueVec[i] <- FALSE
  }
}

par(cex.lab = 1.5)

iterProd <- as.data.frame(feather::read_feather("iterProd.feather"))

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
ATP <- sapply(tails, function(x) if (is.null(x)) {return(NA)} else {unlist(x)[37]})

## Complex I
pdf("../dataVis/complexIATPunivar.pdf")
par(cex.lab = 1.5)
plot(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                  iterProd[ , 1] != 1.0) & univar],
     ylim = c(0, max(ATP, na.rm = T)), ylab = "ATP Concentration",
     xlim = c(0.2, 1.1), xlab = "Proportion of Typical Activity",
     cex = 0.1, main = "Complex I")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
     ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |iterProd[ , 1] != 1.0) & univar],
     col = "springgreen")
abline(a = ATP[which.min(ATP)], b = 0, col = "red")
dev.off()

## Complex III
pdf("../dataVis/complexIIIATPunivar.pdf")
par(cex.lab = 1.5)
plot(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                   iterProd[ , 2] != 1.0) & univar],
     ylim = c(0, max(ATP, na.rm = T)), ylab = "ATP Concentration",
     xlim = c(0.1, 1.1), xlab = "Proportion of Typical Activity", 
     cex = 0.1,
     main = "Complex III")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
     ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 2] != 1.0) & univar],
     col = "springgreen")
abline(a = ATP[which.min(ATP)], b = 0, col = "red")
dev.off()

## Complex IV
pdf("../dataVis/complexIVATPunivar.pdf")
par(cex.lab = 1.5)
plot(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                   iterProd[ , 3] != 1.0) & univar],
     ylim = c(0, max(ATP, na.rm = T)), ylab = "ATP Concentration",
     xlim = c(0.1, 1.1), xlab = "Proportion of Typical Activity", cex = 0.1,
     main = "Complex IV")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
     ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 3] != 1.0) & univar],
     col = "springgreen")
abline(a = ATP[which.min(ATP)], b = 0, col = "red")
dev.off()

## ATP Synthase
pdf("../dataVis/atpsynthaseATPunivar.pdf")
par(cex.lab = 1.5)
plot(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                   iterProd[ , 4] != 1.0) & univar],
     ylim = c(0, max(ATP, na.rm = T)), ylab = "ATP Concentration",
     xlim = c(0.1, 1.1), xlab = "Proportion of Typical Activity",
     cex = 0.1, main = "ATP Synthase")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
     ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 4] != 1.0) & univar],
     col = "springgreen")
abline(a = ATP[which.min(ATP)], b = 0, col = "red")
dev.off()

pdf("../dataVis/univarPlates.pdf")
par(mfrow = c(2, 2), cex.lab = 1.5, cex.main = 1.5)

plot(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                   iterProd[ , 1] != 1.0) & univar],
     ylim = c(0, max(ATP, na.rm = T)), ylab = "[ATP]_c (M)",
     xlim = c(0.2, 1.1), xlab = "Proportion of Typical Activity",
     cex = 0.1, main = "Complex I")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4)-0.0001, c(0.25, 0.5, 0.75, 1)+0.05, 
     ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |iterProd[ , 1] != 1.0) & univar],
     col = "springgreen")
abline(a = ATP[which.min(ATP)], b = 0, col = "red")

plot(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                   iterProd[ , 2] != 1.0) & univar],
     ylim = c(0, max(ATP, na.rm = T)), ylab = "[ATP]_c (M)",
     xlim = c(0.1, 1.1), xlab = "Proportion of Typical Activity", 
     cex = 0.1,
     main = "Complex III")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4)-0.0001, c(0.25, 0.5, 0.75, 1)+0.05, 
     ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 2] != 1.0) & univar],
     col = "springgreen")
abline(a = ATP[which.min(ATP)], b = 0, col = "red")

plot(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                   iterProd[ , 3] != 1.0) & univar],
     ylim = c(0, max(ATP, na.rm = T)), ylab = "[ATP]_c (M)",
     xlim = c(0.1, 1.1), xlab = "Proportion of Typical Activity", 
     cex = 0.1,
     main = "Complex IV")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4)-0.0001, c(0.25, 0.5, 0.75, 1)+0.05, 
     ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 3] != 1.0) & univar],
     col = "springgreen")
abline(a = ATP[which.min(ATP)], b = 0, col = "red")

plot(c(0.25, 0.5, 0.75, 1), ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                   iterProd[ , 4] != 1.0) & univar],
     ylim = c(0, max(ATP, na.rm = T)), ylab = "[ATP]_c (M)",
     xlim = c(0.1, 1.1), xlab = "Proportion of Typical Activity",
     cex = 0.1, main = "ATP Synthase")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(0, 4)-0.0001, c(0.25, 0.5, 0.75, 1)+0.05, 
     ATP[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 4] != 1.0) & univar],
     col = "springgreen")
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
     xlim = c(0.1, 1.1), xlab = "Proportion of Typical Activity",
     cex = 0.1, main = "Complex I")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(149.57, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
     dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 1] != 1.0) & univar],
     col = "springgreen")
abline(a = dPsi[which.min(dPsi)], b = 0, col = "red")
dev.off()

## Ifosfamide

pdf("../dataVis/ifosfamide.pdf")
par(cex.lab = 1.5)
plot(c(0.5, 1), dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
            iterProd[ , 1] == 0.5) & univar],
     ylim = c(155, max(dPsi, na.rm = T)), ylab = "Electrical Potential Gradient",
     xlim = c(0.1, 1.1), xlab = "Proportion of Typical Activity",
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
     xlim = c(0.1, 1.1), xlab = "Proportion of Typical Activity", cex = 0.1,
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
     xlim = c(0.1, 1.1), xlab = "Proportion of Typical Activity", cex = 0.1,
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
     xlim = c(0.1, 1.1), xlab = "Proportion of Typical Activity",
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
     xlim = c(0.1, 1.1), xlab = "Proportion of Typical Activity",
     cex = 0.1, main = "Complex I")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(149.57, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
     dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 1] != 1.0) & univar],
     col = "springgreen")
abline(a = dPsi[which.min(dPsi)], b = 0, col = "red")

plot(c(0.25, 0.5, 0.75, 1), dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                    iterProd[ , 2] != 1.0) & univar],
     ylim = c(150, max(dPsi, na.rm = T)), ylab = "Electrical Potential Gradient",
     xlim = c(0.1, 1.1), xlab = "Proportion of Typical Activity", cex = 0.1,
     main = "Complex III")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(149.57, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
     dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 2] != 1.0) & univar],
     col = "springgreen")
abline(a = dPsi[which.min(dPsi)], b = 0, col = "red")

plot(c(0.25, 0.5, 0.75, 1), dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                    iterProd[ , 3] != 1.0) & univar],
     ylim = c(150, max(dPsi, na.rm = T)), ylab = "Electrical Potential Gradient",
     xlim = c(0.1, 1.1), xlab = "Proportion of Typical Activity", cex = 0.1,
     main = "Complex IV")
rect(c(0.25, 0.5, 0.75, 1)-0.05, rep(149.57, 4), c(0.25, 0.5, 0.75, 1)+0.05, 
     dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 | iterProd[ , 3] != 1.0) & univar],
     col = "springgreen")
abline(a = dPsi[which.min(dPsi)], b = 0, col = "red")

plot(c(0.25, 0.5, 0.75, 1), dPsi[(matrixStats::rowProds(as.matrix(iterProd)) == 1.0 |
                                    iterProd[ , 4] != 1.0) & univar],
     ylim = c(150, max(dPsi, na.rm = T)), ylab = "Electrical Potential Gradient",
     xlim = c(0.1, 1.1), xlab = "Proportion of Typical Activity",
     cex = 0.1, main = "ATP Synthase")
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
LMATP <- lm(allResults$ATP_c ~ allResults$CI + allResults$CIII + allResults$CIV + allResults$ATPSynthase, allResults)
print(summary(LMATP))

LMATPMix <- lm(allResults$ATP_c ~ allResults$CI + allResults$CIII + allResults$CIV + allResults$CI*allResults$CIII*allResults$CIV*allResults$ATPSynthase, allResults)
print(summary(LMATPMix))

library(ggplot2)
container <- allResults
container$CIII <- as.factor(container$CIII)

pdf("../dataVis/atpComplexIIImultivar.pdf")
ggplot(container, aes(x = ATP_c)) +
  geom_histogram(aes(color = CIII, fill = CIII), bins = 24) +
  xlab("[ATP]_c (M)") +
  ylab("Frequency") +
  theme(legend.text = element_text(size = 18), axis.title = element_text(size = 18),
        legend.title = element_text(size = 18), axis.text = element_text(size = 18))
dev.off()

container <- allResults[allResults$CIII == 1, ]
container$CIV <- as.factor(container$CIV)

pdf("atpComplexIVmultivar.pdf")
ggplot(container, aes(x = ATP_c*1000)) +
  geom_histogram(aes(color = CIV, fill = CIV), bins = 24) +
  xlab("ATP Concentration (mmol)") +
  ylab("Frequency") +
  theme(legend.text = element_text(size = 18), axis.title = element_text(size = 18),
        legend.title = element_text(size = 18), axis.text = element_text(size = 18))
dev.off()

# LMATPMix <- lm(allResults$ATP_c ~ allResults$CI + allResults$CIII + allResults$CIV + allResults$ATPSynthase + 
#               allResults$CI*allResults$CIII*allResults$CIV*allResults$ATPSynthase, allResults)
# summary(LMATPMix)

## The dPsi response is well-approximated by a linear surface in CI, CIII, CIV, ATP Synthase, adding mixed terms doesn't
## improve the fit, suggesting that multiple causes interact additively
LMPSI <- lm(allResults$dPsi ~ allResults$CI + allResults$CIII + allResults$CIV, allResults)
print(summary(LMPSI))

LMPSIMix <- lm(allResults$dPsi ~ allResults$CI + allResults$CIII + allResults$CIV, allResults)
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

