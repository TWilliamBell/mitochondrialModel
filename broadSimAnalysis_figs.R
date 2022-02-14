## In order to run this file, you should first run broadSimmTAL.py, writeProduct.py, and then
## collectBroadSim.R

if (!grepl("mitochondrialModel/modelScripts", getwd())) {
  setwd("./modelScripts")
}

drugSimParam <- as.data.frame(feather::read_feather(
  "../results/iterProdBroadSim.feather"))

colnames(drugSimParam) <- c("k", "CI", "CIII", "CIV", "ATPSynthase",
                            "XHle", "glycolysisLevel", "O2Level")

drugSimTails <- read.csv("../results/tailsBroadSimmTAL.csv")
normal <- data.table::fread("../results/resultsmTAL.csv")
normalTail <- tail(normal, 1)
normalTail$V1 <- normalTail$t <- NULL
normalTail <- unlist(normalTail)
drugSimTails$X <- drugSimTails$V1 <- drugSimTails$t <- NULL

mitDis <- data.table::fread("../results/tailsMitDismTALSim.csv")
mitDis$V1 <- mitDis$t <- NULL
colnames(mitDis)[1] <- "k"

iterProd <- feather::read_feather('../results/iterProd.feather')
iterProd <- cbind(1:256, iterProd, rep(1, 256), rep(0, 256), 
                  rep(1, 256))
colnames(iterProd) <- c("k", "CI", "CIII", "CIV", "ATPSynthase",
                        "XHle", "glycolysisLevel", "O2Level")
mitDis <- dplyr::inner_join(iterProd, mitDis)

colnames(drugSimTails)[1] <- "k"
joinedTable <- dplyr::inner_join(drugSimParam, drugSimTails, by = 'k')
mitDis$V1 <- NULL
joinedTable <- rbind(joinedTable, mitDis)
write.csv(joinedTable, "../results/broadSimAllInOneTable.csv")

#joinedTable <- joinedTable[!duplicated(joinedTable[2:8]), ]

## Not especially sensitive to changes due to dichloroacetate alone
O2Only <- joinedTable[joinedTable$CI == 1 & joinedTable$CIII == 1
                       & joinedTable$CIV == 1 
                       & joinedTable$ATPSynthase == 1 
                       & joinedTable$XHle == 1 
                       & joinedTable$glycolysisLevel == 0, ]
## Not very impactful


## Hydrogen leak

XHleOnly <- joinedTable[joinedTable$CI == 1 & joinedTable$CIII == 1
                      & joinedTable$CIV == 1 
                      & joinedTable$ATPSynthase == 1 
                      & joinedTable$glycolysisLevel == 0
                      & joinedTable$O2Level == 1, ]
## Not very impactful on ATPc
uncouplingdPsi <- XHleOnly$dPsi

pdf("../dataVis/uncouplingmTAL.pdf")
par(cex.lab = 1.5, cex.axis = 1.5)
plot(uncouplingdPsi, cex = 0.01, xlim = c(0.5, 4.5), 
     ylim = c(155, 165),
     xaxt = "n", xlab = "Fold Change in Hydrogen Leak", 
     ylab = "Electrical Potential Gradient (mV)")
rect(1:4-0.5, rep(155,4)-0.5, 1:4+0.5, uncouplingdPsi, 
     col = "springgreen")
axis(1, at = 1:4, labels = paste0(XHleOnly$XHle, "x"))
dev.off()

library(ggplot2)
#library(MASS)

pdf("../dataVis/dPsileakOXPHOSO2mTALmultivar.pdf")

container <- joinedTable
container$O2Level <- as.factor(joinedTable$O2Level)
p <- ggplot(container, aes(x = dPsi)) +
  geom_histogram(aes(fill = O2Level)) +
  xlab("dPsi (mV)") +
  ylab("Frequency") +
  theme(axis.title = element_text(size = 18), legend.title = element_text(size = 18),
        legend.text = element_text(size = 18), axis.text = element_text(size = 18))
p$labels$fill <- "Relative \nComplex IV \nActivity"
p

dev.off()

## Figure 3.16
pdf("../dataVis/atpleakOXPHOSO2mTALmultivar.pdf")
container <- joinedTable
container <- container[container$glycolysisLevel == 0., ]
container$O2Level <- as.factor(container$O2Level)
p <- ggplot(container, aes(x = ATP_c*1000)) +
  geom_histogram(aes(fill = O2Level)) +
  xlab("Cytosolic ATP (mM)") +
  ylab("Frequency") +
  xlim(c(0, 2.5))
p$labels$fill <- "Relative \nOxygen \nTension"
p +
  theme(axis.title = element_text(size = 18), legend.title = element_text(size = 18),
        legend.text = element_text(size = 18), axis.text = element_text(size = 18))
dev.off()

container <- joinedTable[joinedTable$glycolysisLevel == 0 & 
                           joinedTable$O2Level == 1 &
                           joinedTable$XHle == 1,]
container$CIV <- as.factor(container$CIV)
ggplot(container, aes(x = dPsi)) +
  geom_histogram(aes(color = CIV, fill = CIV)) +
  xlab("dPsi (mV)") +
  ylab("Frequency")

## Figure 3.15
pdf("../dataVis/atpCIVmTAL.pdf")
container <- mitDis
container$CIV <- as.factor(container$CIV)
par(cex.lab = 1.5, cex.axis = 1.5)
p <- ggplot(container, aes(x = ATP_c*1000)) +
  geom_histogram(aes(fill = CIV), bins = 20) +
  xlab("Cytosolic ATP (mM)") +
  ylab("Frequency") + 
  theme(axis.title = element_text(size = 18), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18), 
        axis.text = element_text(size = 18))
p$labels$fill <- "Relative\nComplex IV\nActivity"
p
dev.off()

pdf("../dataVis/atpHist.pdf")
par(cex.lab = 1.5, cex.axis = 1.5)
hist(joinedTable$ATP_c, main = NULL, xlab = "[ATP]_c (M)")
dev.off()

joinedTable$glycolysisLevel <- joinedTable$glycolysisLevel/max(joinedTable$glycolysisLevel)

## Table 3.7 uses these values
LMdPsi <- lm(dPsi ~ CI + CIII + CIV + XHle + O2Level + glycolysisLevel, 
             joinedTable)
summary(LMdPsi)

LMATP <- lm(ATP_c ~ CI + CIII + CIV + XHle + O2Level + glycolysisLevel 
            + CIV*XHle*O2Level,
            joinedTable)
summary(LMATP)

## OXPHOS univariate
dPsiMin <- min(mitDis$dPsi)

vals <- (joinedTable$CI == 1) + (joinedTable$CIII == 1) + 
  (joinedTable$CIV == 1) + (joinedTable$ATPSynthase == 1)
vals <- vals == 3

oxphos <- joinedTable[joinedTable$XHle == 1 &
                        joinedTable$glycolysisLevel == 0 &
                        joinedTable$O2Level == 1 & vals, ]

#noms <- c("CI", "CIII", "CIV", "ATPSynthase")

pdf("../dataVis/dPsiUnivarmTAL.pdf")

par(cex.axis = 1.5)
par(mfrow = c(2, 2), mar = c(2, 2, 1, 1))
for (i in 0:3) {
  relevantBit <- oxphos[3*i+1:3, ]
  dPsi <- c(relevantBit$dPsi, unname(normalTail[2]))
  #ATP <- c(relevantBit$ATP_c, unname(normalTail[35]))
  #pdf(paste0("./dataVis/univardPsimTAL", noms[i+1], ".pdf"))
  plot(dPsi, cex = 0.01, xlim = c(0.5, 4.5), ylim = c(155, 170),
       xlab = '', xaxt = "n")
  rect(1:4-0.25, rep(155-0.6, 4), 1:4+0.25, dPsi, col = "springgreen")
  axis(side = 1, at = 1:4, labels = c(0.25, 0.5, 0.75, 1))
  abline(a = dPsiMin, b = 0, col = "red")
  #dev.off()
  
  #pdf(paste0("./results/univarATPmTAL", noms[i+1], ".pdf"))
  #plot(ATP, cex = 0.01, xlim = c(0.5, 4.5), ylim = c(0, 0.003))
  #rect(1:4-0.5, rep(0-0.00013, 4), 1:4+0.5, ATP, col = "springgreen")
  #dev.off()
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

## Mitochondrial disease multivariate
container <- mitDis
container$CIII <- as.factor(container$CIII)
container$CIV <- as.factor(container$CIV)

pdf("../dataVis/dPsiMitDismTAL.pdf")
ggplot(container, aes(x = dPsi)) +
  geom_histogram(aes(color = CIII, fill = CIII)) +
  xlab("dPsi (mV)") +
  ylab("Frequency") + 
  theme(axis.title = element_text(size = 18), legend.title = element_text(size = 18),
        legend.text = element_text(size = 18), axis.text = element_text(size = 18))
dev.off()

# ggplot(container, aes(x = dPsi)) +
#   geom_histogram(aes(color = CIV, fill = CIV)) +
#   xlab("dPsi (mV)") +
#   ylab("Frequency")

## Worst cases for ATP
#joinedTable[joinedTable$ATP_c < 0.002, ]