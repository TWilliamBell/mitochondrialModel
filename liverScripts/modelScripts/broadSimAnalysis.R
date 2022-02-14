if (!grepl("liverScripts/modelScripts", getwd())) {
  setwd("./liverScripts/modelScripts")
}

tails <- read.csv("../results/tailsBroadSim.csv")
rownames(tails) <- NULL

iterprod <- read.csv("../results/iterprod.csv")

iterprod$X <- NULL
tails$X <- tails$V1 <- tails$t <- NULL

broadSim <- dplyr::inner_join(iterprod, tails, by = "nums")

oxphos <- broadSim[broadSim$Glycolysis == 0.00023 & broadSim$PO2 == 1.0, ]

LMoxphosATP <- lm(oxphos$ATP_c ~ oxphos$CI + oxphos$CIII + oxphos$CIV + 
                    oxphos$ATP)# + oxphos$CI*oxphos$CIII*oxphos$CIV*oxphos$ATP)
summary(LMoxphosATP)

LMoxphosdPsi <- lm(oxphos$dPsi ~ oxphos$CI + oxphos$CIII + oxphos$CIV + 
                    oxphos$ATP)# + oxphos$CI*oxphos$CIII*oxphos$CIV*oxphos$ATP)
summary(LMoxphosdPsi)

CI <- oxphos[oxphos$CIII == 1.0 & oxphos$CIV == 1.0 & oxphos$ATP == 1.0, ]
CIII <- oxphos[oxphos$CI == 1.0 & oxphos$CIV == 1.0 & oxphos$ATP == 1.0, ]
CIV <- oxphos[oxphos$CI == 1.0 & oxphos$CIII == 1.0 & oxphos$ATP == 1.0, ]
ATP <- oxphos[oxphos$CIII == 1.0 & oxphos$CIV == 1.0 & oxphos$CI == 1.0, ]

barplot(CI$dPsi, xlab = "Normalized Complex I Activity", 
        names.arg = c(0.25, 0.5, 0.75, 1.0), col = "thistle4", xpd = F,
        ylim = c(140, 170))
barplot(CIII$dPsi, xlab = "Normalized Complex III Activity",
        names.arg = c(0.25, 0.5, 0.75, 1.0), col = "thistle4", xpd = F,
        ylim = c(140, 170))
barplot(CIV$dPsi, xlab = "Normalized Complex IV Activity",
        names.arg = c(0.25, 0.5, 0.75, 1.0), col = "thistle4", xpd = F,
        ylim = c(140, 170))
barplot(ATP$dPsi, xlab = "Normalized ATP Synthase Activity",
        names.arg = c(0.25, 0.5, 0.75, 1.0), col = "thistle4", xpd = F,
        ylim = c(140, 170))

barplot(CI$ATP_c, xlab = "Normalized Complex I Activity", 
        names.arg = c(0.25, 0.5, 0.75, 1.0), col = "thistle4")
barplot(CIII$ATP_c, xlab = "Normalized Complex III Activity",
        names.arg = c(0.25, 0.5, 0.75, 1.0), col = "thistle4")
barplot(CIV$ATP_c, xlab = "Normalized Complex IV Activity",
        names.arg = c(0.25, 0.5, 0.75, 1.0), col = "thistle4")
barplot(ATP$ATP_c, xlab = "Normalized ATP Synthase Activity",
        names.arg = c(0.25, 0.5, 0.75, 1.0), col = "thistle4")

library(ggplot2)

container <- oxphos
container$CIII <- as.factor(container$CIII)
ggplot(container, aes(x = dPsi)) +
  geom_histogram(aes(color = CIII, fill = CIII)) +
  xlab("dPsi (mV)") +
  ylab("Frequency")

container <- oxphos
container$CIII <- as.factor(container$CIII)
ggplot(container, aes(x = 1000*ATP_c)) +
  geom_histogram(aes(color = CIII, fill = CIII)) +
  xlab("ATP (mM)") +
  ylab("Frequency")

## Run through it again with P_O2 and glycolysis

container <- broadSim
container$PO2 <- as.factor(container$PO2)
ggplot(container, aes(x = dPsi)) +
  geom_histogram(aes(color = PO2, fill = PO2)) +
  xlab("dPsi (mV)") +
  ylab("Frequency")

container <- broadSim
container$PO2 <- as.factor(container$PO2)
ggplot(container, aes(x = ATP_c)) +
  geom_histogram(aes(color = PO2, fill = PO2)) +
  xlab("ATP_c (M)") +
  ylab("Frequency")

container <- broadSim
container$Glycolysis <- as.factor(container$Glycolysis)
ggplot(container, aes(x = 1000*ATP_c)) +
  geom_histogram(aes(color = Glycolysis, fill = Glycolysis)) +
  xlab("ATP (mM)") +
  ylab("Frequency")


PO2 <- broadSim[broadSim$CI == 1.0 & broadSim$CIII == 1.0 & broadSim$CIV == 1.0 & 
                  broadSim$ATP == 1.0 & broadSim$Glycolysis == 0.00023, ]

Glycolysis <- broadSim[broadSim$CI == 1.0 & broadSim$CIII == 1.0 & 
                         broadSim$CIV == 1.0 & broadSim$ATP == 1.0 & 
                         broadSim$PO2 == 1.0, ]

broadSim$Glycolysis <-  round(broadSim$Glycolysis/0.00023, 1)
LMbroadATP <- lm(broadSim$ATP_c ~ broadSim$CI + broadSim$CIII + broadSim$CIV + 
                    broadSim$ATP + broadSim$Glycolysis + broadSim$PO2)
summary(LMbroadATP)

LMbroaddPsi <- lm(broadSim$dPsi ~ broadSim$CI + broadSim$CIII + broadSim$CIV + 
                     broadSim$ATP + broadSim$Glycolysis + broadSim$PO2)
summary(LMbroaddPsi)

