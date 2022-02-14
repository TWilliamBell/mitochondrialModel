if (!grepl("liverScripts/modelscripts", getwd())) {
  setwd("./liverScripts/modelScripts")
}

tails <- read.csv("../results/tailsDiabetes.csv")
rownames(tails) <- NULL

iterprod <- read.csv("../results/iterprodDiabetes.csv")

diabetes <- dplyr::inner_join(iterprod, tails, by = "nums")
diabetes$X.x <- diabetes$X.y <- diabetes$V1 <- diabetes$t <- NULL

Leak <- diabetes[diabetes$CI == 1.0 & diabetes$CIII == 1.0 & diabetes$CIV == 1.0 &
                   diabetes$ATP == 1.0 & diabetes$PO2 == 1.0, ]

barplot(Leak$dPsi, 
        names.arg = c(1.0, 1.15, 1.5, 5.0, 10.0), col = "thistle4",
        xpd = F,
        ylim = c(150, 165))
barplot(Leak$ATP_c, xpd = F, ylim = c(0.006, 0.007),
        names.arg = c(1.0, 1.15, 1.5, 5.0, 10.0), col = "thistle4")
# barplot(Leak$Cred_i/0.0058387, names.arg = c(1.0, 1.15, 1.5, 5.0, 10.0), col = "thistle4",
#         xpd = F)

PO2 <- diabetes[diabetes$CI == 1.0 & diabetes$CIII == 1.0 & diabetes$CIV == 1.0 &
                  diabetes$ATP == 1.0 & diabetes$Leak == 1.0, ]

barplot(PO2$dPsi, 
        names.arg = c(0.1, 0.5, 1.0), col = "thistle4",
        xpd = F,
        ylim = c(150, 165))
barplot(PO2$ATP_c, xpd = F, ylim = c(0.006, 0.007),
        names.arg = c(0.1, 0.5, 1.0), col = "thistle4")

library(ggplot2)

container <- diabetes
container$PO2 <- as.factor(container$PO2)
ggplot(container, aes(x = dPsi)) +
  geom_histogram(aes(color = PO2, fill = PO2)) +
  xlab("dPsi (mV)") +
  ylab("Frequency")

pdf("../dataVis/atpLiverDiabetes.pdf")
container <- diabetes
container$Leak <- as.factor(container$Leak)
ggplot(container, aes(x = 1000*ATP_c)) +
  geom_histogram(aes(color = Leak, fill = Leak)) +
  xlab("ATP (mM)") +
  ylab("Frequency") + 
  theme(axis.title = element_text(size = 18), legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))
dev.off()

LeakPO2 <- diabetes[diabetes$CI == 1.0 & diabetes$CIII == 1.0 &
                      diabetes$CIV == 1.0 & diabetes$ATP == 1.0, ]

container <- LeakPO2
container$Leak <- as.factor(container$Leak)
ggplot(container, aes(x = dPsi)) +
  geom_histogram(aes(color = Leak, fill = Leak)) +
  xlab("dPsi (mV)") +
  ylab("Frequency")

container <- LeakPO2
container$Leak <- as.factor(container$Leak)
ggplot(container, aes(x = 1000*ATP_c)) +
  geom_histogram(aes(color = Leak, fill = Leak)) +
  xlab("ATP (mM)") +
  ylab("Frequency")

LeakPO2$Leak <- as.factor(LeakPO2$Leak)
ggplot(LeakPO2, aes(x = Cred_i/0.0058387)) +
  geom_histogram(aes(color = Leak, fill = Leak)) +
  xlab("Cytochrome C Redox State") +
  ylab("Frequency")

oxphosOnly <- diabetes[diabetes$Leak == 1 & diabetes$PO2 == 1, ]
oxphosOnly$CIV <- as.factor(oxphosOnly$CIV)
ggplot(oxphosOnly, aes(x = ATP_c*1000)) +
  geom_histogram(aes(color = CIV, fill = CIV)) +
  xlab("ATP (mM)") +
  ylab("Frequency")

ggplot(oxphosOnly, aes(x = dPsi)) +
  geom_histogram(aes(color = CIV, fill = CIV)) +
  xlab("dPsi (mV)") +
  ylab("Frequency")

ggplot(diabetes, aes(x = QH2_x/0.001298)) +
  geom_histogram(aes(fill = Leak)) +
  xlab("Coenzyme Q Reduction State")

pdf("../dataVis/cytCdiabetesLiver.pdf")
diabetes$CIV <- as.factor(diabetes$CIV)
ggplot(diabetes, aes(x = Cred_i/0.0058387)) +
  geom_histogram(aes(colour = CIV, fill = CIV)) +
  xlab("Cytochrome C Redox State") +
  ylab("Frequency") +
  geom_vline(xintercept = diabetes[423, ]$Cred_i/0.0058387, col = "red") +
  xlim(c(0, 1)) + 
  theme(axis.title = element_text(size = 18), legend.title = element_text(size = 18),
        legend.text = element_text(size = 18), axis.text = element_text(size = 18))
dev.off()

ggplot(diabetes, aes(x = NADH_x/9.4e-3)) +
  geom_histogram(aes(colour = Leak, fill = Leak)) +
  xlab("NADH/NAD+ Reduction State")

