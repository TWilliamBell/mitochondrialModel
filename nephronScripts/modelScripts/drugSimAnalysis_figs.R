if (!grepl("nephronScripts/modelScripts", getwd())) {
  setwd("./nephronScripts/modelScripts")
}

## Before running this script, run drugSim.py & writeProduct.py and then collectDrugSim.R

drugSimParam <- data.table::fread("../results/iterProdDrugSim.csv")
drugSimParam$V1 <- NULL

colnames(drugSimParam) <- c("k", "CI", "CIII", "CIV", "ATPSynthase",
                            "XHle", "PDHEq")

drugSimTails <- read.csv("../results/tailsDrugSim.csv")
normal <- data.table::fread("../results/resultsATP.csv")
normalTail <- tail(normal, 1)
normalTail$V1 <- normalTail$t <- NULL
normalTail <- unlist(normalTail)
drugSimTails$X <- drugSimTails$V1 <- drugSimTails$t <- NULL


colnames(drugSimTails)[1] <- "k"
joinedTable <- dplyr::inner_join(drugSimParam, drugSimTails, by = 'k')
drugSimParam <- drugSimParam[!(drugSimParam$PDHEq == 1
              & drugSimParam$XHle == 1), ]

#joinedTable <- cbind(drugSimParam, drugSimTails)
## Not especially sensitive to changes due to dichloroacetate alone
DCAOnly <- joinedTable[joinedTable$CI == 1 & joinedTable$CIII == 1 & 
                         joinedTable$CIV == 1 & joinedTable$ATPSynthase == 1 
                       & joinedTable$XHle == 1, ]
DCARatio <- apply(as.matrix(DCAOnly)[ , -(1:7)], 1, 
                  function(x) x/normalTail)

## Most notable change that there are higher concentrations of most TCA
## cycle intermediates, what this says is that DCA is not super
## impactful on its own.

## Ifosfamide alone covered in atpResponseMitDis.R.  Now we consider
## DCA with Ifosfamide.

# dcadPsi <- c(joinedTable$dPsi[joinedTable$CI == 0.5 & 
#                   joinedTable$CIII == 1 & 
#                   joinedTable$CIV == 1 & 
#                   joinedTable$ATPSynthase == 1 & 
#                   joinedTable$XHle == 1][1],
#              159.6287,
#              159.8163, 
#              160.3137) ## Numbers taken from mit dis, DCA

dcaATP <- c(joinedTable$ATP_c[joinedTable$CI == 0.5 & 
                               joinedTable$CIII == 1 & 
                               joinedTable$CIV == 1 & 
                               joinedTable$ATPSynthase == 1 & 
                               joinedTable$XHle == 1][1],
            0.002318198,
            0.002330648, 
            0.002359787) ## See above

# pdf("../dataVis/ifosfamidedPsi.pdf")
# par(cex.lab = 1.5)
# plot(dcadPsi, ylim = c(155, max(dcadPsi)), xlim = c(0, 5),
#      xaxt = "n", xlab = "Conditions", cex = 0.1, ylab = "dPsi (mV)")
# rect(1:4-0.5, rep(154.78, 4), 1:4+0.5, dcadPsi, col = "springgreen")
# axis(side = 1, at = 1:4, labels = c("DCA+/Ifos+", "DCA+/Ifos-", 
#                                     "DCA-/Ifos+", "DCA-/Ifos-"))
# dev.off()

pdf("../dataVis/ifosfamideATP.pdf", width = 12)
par(cex.lab = 2.5, cex.axis = 1.75, mar = c(5.1, 5.1, 4.1, 2.1))
plot(dcaATP*1000, ylim = c(0.0001, max(dcaATP))*1000, xlim = c(0, 5),
     xaxt = "n", xlab = "Conditions", cex = 0.1, ylab = "Cytosolic ATP (mM)")
rect(1:4-0.5, rep(0, 4), 1:4+0.5, dcaATP*1000, col = "springgreen")
axis(side = 1, at = 1:4, labels = c("DCA+/Ifos+", "DCA+/Ifos-", 
                "DCA-/Ifos+", "DCA-/Ifos-"))
dev.off()

## Hydrogen Leak
SalicylateOnly <- joinedTable[joinedTable$CI == 1 & joinedTable$CIII == 1 & 
                         joinedTable$CIV == 1 & joinedTable$ATPSynthase == 1 
                       & joinedTable$PDHEq == 1, ]
dPsiSal <- c(160.3137, SalicylateOnly$dPsi)

pdf("../dataVis/hleakRange.pdf")
par(cex.lab = 1.5)
plot(dPsiSal, ylim = c(155, max(dPsiSal)), xlim = c(0.5, 6.5),
     xaxt = "n", xlab = "Conditions", cex = 0.1, ylab = "dPsi (mV)")
rect(1:6-0.5, rep(154.78, 5), 1:6+0.5, dPsiSal, col = "springgreen")
axis(side = 1, at = 1:6, labels = c("Norm H-Leak", "1.15x", "1.5x",
                                    '2x', "5x", "10x"))
dev.off()

ATPSal <- c(0.0023, SalicylateOnly$ATP_c)
pdf("../dataVis/hleakRangeATP.pdf")
par(cex.lab = 1.5)
plot(ATPSal, ylim = c(0, max(ATPSal)), xlim = c(0.5, 6.5),
     xaxt = "n", xlab = "Conditions", cex = 0.1, ylab = "ATP (mM)")
rect(1:6-0.5, rep(0, 5), 1:6+0.5, ATPSal, col = "springgreen")
axis(side = 1, at = 1:6, labels = c("Norm H-Leak", "1.15x", "1.5x",
                                    '2x', "5x", "10x"))
dev.off()

## Hydrogen Leak with DCA
SalicylateDCA <- joinedTable[joinedTable$CI == 1 & joinedTable$CIII == 1 & 
                                joinedTable$CIV == 1 & joinedTable$ATPSynthase == 1 
                              & joinedTable$PDHEq == 2.0, ]
dPsiSal <- c(160.3137, SalicylateDCA$dPsi)

# plot(dPsiSal, ylim = c(155, max(dPsiSal)), xlim = c(0.5, 7.5),
#      xaxt = "n", xlab = "Conditions", cex = 0.1, ylab = "dPsi (mV)")
# rect(1:7-0.5, rep(154.78, 5), 1:7+0.5, dPsiSal, col = "springgreen")
# axis(side = 1, at = 1:6, labels = c("Norm H-Leak", "1.15x", "1.5x",
#                                     '2x', "5x", "10x"))
## Does not appear significantly different from DCA- case

## Oxphos disorder plus Hleak

mitDis <- read.csv("../results/tailsMitDisPT.csv")
mitDis <- mitDis[order(mitDis$nums), ]
mitDisIter <- data.table::fread("../results/iterProd.csv")
mitDisIter$V1 <- NULL
mitDisIter <- mitDisIter[-1, ]
mitDis <- cbind(mitDisIter, mitDis)
mitDisPlusLeak <- joinedTable[joinedTable$PDHEq == 1 , ]
mitDisPlusLeak$k <- mitDisPlusLeak$PDHEq <- NULL
colnames(mitDis)[1:4] <- c("CI", "CIII", "CIV", "ATPSynthase")
mitDis$XHle <- rep(1, 256)
mitDis1 <- mitDis
mitDis1$X <- mitDis1$XHle
mitDis1$XHle <- mitDis1$nums <- mitDis1$V1 <- mitDis1$t <- NULL
colnames(mitDis1)[5] <- "XHle"
#mitDis <- mitDis[ , colnames(mitDisPlusLeak)]
mitDisPlusLeak <- rbind(mitDis1[complete.cases(mitDis1), ], 
                        mitDisPlusLeak)
mitDis1$PDHEq <- rep(1, 256)
rm(mitDis)

library(ggplot2)

LMdPsi <- lm(dPsi ~ CI + CIII + CIV + ATPSynthase + XHle +
               CI*CIII*CIV*ATPSynthase*XHle, mitDisPlusLeak)
summary(LMdPsi)

LMATP <- lm(ATP_c ~ CI + CIII + CIV + ATPSynthase + XHle +
              CI*CIII*CIV*ATPSynthase*XHle, mitDisPlusLeak)
summary(LMATP)

pdf("../dataVis/atpOXPHOSmultivarCIII.pdf")
container <- mitDis1
container$CIII <- as.factor(mitDis1$CIII)
ggplot(container, aes(x = ATP_c)) +
  geom_histogram(aes(color = CIII, fill = CIII)) +
  xlab("ATP_c (M)") +
  ylab("Frequency") +
  theme(legend.text = element_text(size = 18), axis.title = element_text(size = 18),
        legend.title = element_text(size = 18), axis.text = element_text(size = 18))
dev.off()

## Figure 3.13
pdf("../dataVis/dPsileakOXPHOSCIIImultivar.pdf")
container <- mitDisPlusLeak
container$CIII <- as.factor(mitDisPlusLeak$CIII)
p <- ggplot(container, aes(x = dPsi)) +
  geom_histogram(aes(fill = CIII)) +
  xlab("dPsi (mV)") +
  ylab("Frequency") +
  theme(legend.text = element_text(size = 18), axis.title = element_text(size = 18),
        legend.title = element_text(size = 18), axis.text = element_text(size = 18))
p$labels$fill <- "Relative\nComplex III\nActivity"
p
dev.off()

## Figure 3.12
pdf("../dataVis/atpleakOXPHOSCIIImultivar.pdf")
container <- mitDisPlusLeak
container <- container[!(container$CIII == 0.25 & container$ATP_c > 1.8e-3), ]
container$CIII <- as.factor(container$CIII)
p <- ggplot(container, aes(x = 1000*ATP_c)) +
  geom_histogram(aes(fill = CIII)) +
  xlab("Cytosolic ATP (mM)") +
  ylab("Frequency") +
  theme(legend.text = element_text(size = 18), axis.title = element_text(size = 18),
        legend.title = element_text(size = 18), axis.text = element_text(size = 18)) +
  xlim(c(0, 2.5))
p$labels$fill <- "Relative\nComplex III\nActivity"
p
dev.off()

# container <- mitDisPlusLeak[mitDisPlusLeak$CIII == 1.0, ]
# container$CIV <- as.factor(container$CIV)
# ggplot(container, aes(x = dPsi)) +
#   geom_histogram(aes(color = CIV, fill = CIV)) +
#   xlab("dPsi (mV)") +
#   ylab("Frequency")

# container <- mitDisPlusLeak[mitDisPlusLeak$CIII == 1.0
#                             & mitDisPlusLeak$CIV == 1.0, ]
# container$XHle <- as.factor(container$XHle)
# ggplot(container, aes(x = dPsi)) +
#   geom_histogram(aes(color = XHle, fill = XHle)) +
#   xlab("dPsi (mV)") +
#   ylab("Frequency")

## Oxphos disorder, Hleak and DCA

jtBigger <- joinedTable
jtBigger$k <- NULL
mitDis1 <- mitDis1[ , colnames(jtBigger)]
jtBigger <- rbind(jtBigger, mitDis1[complete.cases(mitDis1), ])

LMdPsi <- lm(dPsi ~ CI + CIII + CIV + XHle + PDHEq, jtBigger)
summary(LMdPsi)
LMATP <- lm(ATP_c ~ CI + CIII + CIV + XHle + PDHEq, jtBigger)
summary(LMATP)
#LMATP <- lm(ATP_c ~ CI + CIII + CIV + XHle + PDHEq + CI*CIII*CIV*XHle*PDHEq,
#            jtBigger)
#summary(LMATP)

pdf("../dataVis/dPsileakOXPHOSCIIIDCAmultivar.pdf")
container <- jtBigger
container$CIII <- as.factor(jtBigger$CIII)
ggplot(container, aes(x = dPsi)) +
  geom_histogram(aes(color = CIII, fill = CIII)) +
  xlab("dPsi (mV)") +
  ylab("Frequency") +
  theme(legend.text = element_text(size = 18), axis.title = element_text(size = 18),
        legend.title = element_text(size = 18), axis.text = element_text(size = 18))
dev.off()

pdf("../dataVis/atpleakOXPHOSCIIIDCAmultivar.pdf")
container <- jtBigger
container <- container[!(container$CIII == 0.25 & container$ATP_c > 1.8e-3), ]
container$CIII <- as.factor(container$CIII)
p <- ggplot(container, aes(x = ATP_c*1000)) +
  geom_histogram(aes(fill = CIII)) +
  xlab("Cytosolic ATP (mM)") +
  ylab("Frequency") +
  theme(legend.text = element_text(size = 18), axis.title = element_text(size = 18),
        legend.title = element_text(size = 18), axis.text = element_text(size = 18))
p$labels$fill <- "Relative\nComplex III\nActivity"
p
dev.off()

# container <- jtBigger
# container$PDHEq <- as.factor(jtBigger$PDHEq)
# ggplot(container, aes(x = dPsi)) +
#   geom_histogram(aes(color = PDHEq, fill = PDHEq)) +
#   xlab("dPsi (mV)") +
#   ylab("Frequency")

