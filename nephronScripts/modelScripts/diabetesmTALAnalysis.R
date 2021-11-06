if (!grepl("mitochondrialModel/modelScripts", getwd())) {
  setwd("./modelScripts")
}

tails <- read.csv("../results/tailsDiabetesComplexmTAL.csv")
rownames(tails) <- NULL

iterprod <- feather::read_feather("../results/iterProdDiabetes.feather")
colnames(iterprod) <- c("nums", "CI", "CIII", "CIV", "ATP", "Leak", "PO2")

diabetes <- dplyr::inner_join(iterprod, tails, by = "nums")
diabetes$V1 <- diabetes$t <- diabetes$X <- NULL

Leak <- diabetes[diabetes$CI == 1.0 & diabetes$CIII == 1.0 & diabetes$CIV == 1.0 &
                   diabetes$ATP == 1.0 & diabetes$PO2 == 1.0, ]

barplot(Leak$dPsi, 
        names.arg = c(1.0, 1.15, 1.5, 5.0, 10.0), col = "thistle4",
        xpd = F,
        ylim = c(140, 165))
barplot(Leak$ATP_c, xpd = F, ylim = c(0.00, 0.003),
        names.arg = c(1.0, 1.15, 1.5, 5.0, 10.0), col = "thistle4")
barplot(Leak$QH2_x/0.00178, names.arg = c(1.0, 1.15, 1.5, 5.0, 10.0), col = "thistle4",
        xpd = F)

PO2 <- diabetes[diabetes$CI == 1.0 & diabetes$CIII == 1.0 & diabetes$CIV == 1.0 &
                  diabetes$ATP == 1.0 & diabetes$Leak == 1.0, ]

barplot(PO2$dPsi, 
        names.arg = c(0.1, 0.5, 1.0), col = "thistle4",
        xpd = F,
        ylim = c(140, 165))
barplot(PO2$ATP_c, xpd = F, ylim = c(0.00, 0.003),
        names.arg = c(0.1, 0.5, 1.0), col = "thistle4")

library(ggplot2)

container <- diabetes
container$PO2 <- as.factor(container$PO2)
ggplot(container, aes(x = dPsi)) +
  geom_histogram(aes(color = PO2, fill = PO2)) +
  xlab("dPsi (mV)") +
  ylab("Frequency")

## ATP across all diabetic cases can be depleted
pdf("../dataVis/diabetesmTALATPPO2.pdf")
container <- diabetes
container$PO2 <- as.factor(container$PO2)
ggplot(container, aes(x = 1000*ATP_c)) +
  geom_histogram(aes(color = PO2, fill = PO2)) +
  xlab("ATP (mM)") +
  ylab("Frequency") + 
  theme(axis.title = element_text(size = 18), legend.title = element_text(size = 18),
        legend.text = element_text(size = 18), axis.text = element_text(size = 18))
dev.off()

LeakPO2 <- diabetes[diabetes$CI == 1.0 & diabetes$CIII == 1.0 &
                      diabetes$CIV == 1.0 & diabetes$ATP == 1.0, ]

container <- LeakPO2
container$PO2 <- as.factor(container$PO2)
ggplot(container, aes(x = dPsi)) +
  geom_histogram(aes(color = PO2, fill = PO2)) +
  xlab("dPsi (mV)") +
  ylab("Frequency")

container <- LeakPO2
container$PO2 <- as.factor(container$PO2)
ggplot(container, aes(x = 1000*ATP_c)) +
  geom_histogram(aes(color = PO2, fill = PO2)) +
  xlab("ATP (mM)") +
  ylab("Frequency")

LeakPO2$PO2 <- as.factor(LeakPO2$PO2)
ggplot(LeakPO2, aes(x = Cred_i/2.148e-3)) +
  geom_histogram(aes(color = PO2, fill = PO2)) +
  xlab("Cytochrome C Redox State") +
  ylab("Frequency") +
  xlim(c(0, 1)) +
  geom_vline(xintercept = diabetes[348, ]$Cred_i/2.148e-3, col = "red")

LeakPO2$PO2 <- as.factor(LeakPO2$PO2)
ggplot(LeakPO2, aes(x = QH2_x/0.00178)) +
  geom_histogram(aes(color = PO2, fill = PO2)) +
  xlab("Cytochrome C Redox State") +
  ylab("Frequency") +
  xlim(c(0, 1)) +
  geom_vline(xintercept = diabetes[348, ]$QH2_x/0.00178, col = "red")

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

## CoQ pool left in a highly reduced state under just OXPHOS dysfunction
pdf("../dataVis/oxphosmTALCIVcoq.pdf")
oxphosOnly$CIV <- as.factor(oxphosOnly$CIV)
ggplot(oxphosOnly, aes(x = QH2_x/0.00178)) +
  geom_histogram(aes(color = CIV, fill = CIV)) +
  xlab("Coenzyme Q Redox State") +
  xlim(c(0, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = oxphosOnly[24, ]$QH2_x/0.00178, col = "red") + 
  theme(axis.title = element_text(size = 18), legend.title = element_text(size = 18),
        legend.text = element_text(size = 18), axis.text = element_text(size = 18))
dev.off()

## CoQ pool left in a highly reduced state under combined dysfunction
pdf("../dataVis/diabetesmTALCOQ.pdf")
diabetes$PO2 <- as.factor(diabetes$PO2)
diabetes$CIV <- as.factor(diabetes$CIV)
ggplot(diabetes, aes(x = QH2_x/0.00178)) +
  geom_histogram(aes(color = PO2, fill = PO2)) +
  xlab("Coenzyme Q Redox State") +
  xlim(c(0, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = diabetes[348, ]$QH2_x/0.00178, col = "red") + 
  theme(axis.title = element_text(size = 18), legend.title = element_text(size = 18),
        legend.text = element_text(size = 18), axis.text = element_text(size = 18))
dev.off()

pdf("../dataVis/diabetesmTALCytC.pdf")
ggplot(diabetes, aes(x = Cred_i/2.148e-3)) +
  geom_histogram(aes(color = PO2, fill = PO2)) +
  xlab("Cytochrome C Redox State") +
  ylab("Frequency") +
  xlim(c(0, 1)) +
  geom_vline(xintercept = diabetes[348, ]$Cred_i/2.148e-3, col = "red") + 
  theme(axis.title = element_text(size = 18), legend.title = element_text(size = 18),
        legend.text = element_text(size = 18), axis.text = element_text(size = 18))
dev.off()

ggplot(diabetes, aes(x = NADH_x/0.824e-4)) +
  geom_histogram() +
  xlab("NADH/NAD+ Redox State") +
  ylab("Frequency") +
  xlim(c(0, 1)) +
  geom_vline(xintercept = diabetes[348, ]$NADH_x/0.824e-4, col = "red")



uncoupling <- paste0("../results/", dir("../results")[grepl("MiceUncouplingmTAL", 
                                                      dir("../results"))])
newCases <- list()

for (i in seq_along(uncoupling)) {
  newCases[[i]] <- data.table::fread(uncoupling[i])
}

#atp <- function(x) tail(x$ATP_c, 1)
#newCasesATP <- sapply(newCases, atp)
#hist(newCasesATP/0.00258)

cyt <- function(x) tail(x$Cred_i, 1)
newCasesCytC <- sapply(newCases, cyt)
hist(newCasesCytC/2.148e-3, xlim = c(0, 1))

coq <- function(x) tail(x$QH2_x, 1)
newCasesCOQ <- sapply(newCases, coq)
hist(newCasesCOQ/0.00178)

nadh <- function(x) tail(x$NADH_x, 1)
newCasesNADH <- sapply(newCases, nadh)
hist(newCasesNADH/0.824e-4)
