## In order to run this file, you should first run (in kidney and liver directories) diabetesComplex.py, 
## ("") diabetesCollectVals.R, and ("") writeProduct.py

# if (!grepl("mitochondrialModel/nephronScripts/modelScripts", getwd())) {
#   setwd("./nephronScripts/modelScripts")
# }
setwd("modelScripts")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

tails <- read.csv("../results/tailsDiabetesComplex.csv")
rownames(tails) <- NULL

iterprod <- feather::read_feather("../results/iterProdDiabetes.feather")
colnames(iterprod) <- c("nums", "CI", "CIII", "CIV", "ATP", "Leak", "PO2")

diabetesPT <- dplyr::inner_join(iterprod, tails, by = "nums")
diabetesPT$V1 <- diabetesPT$t <- diabetesPT$X <- NULL

tails <- read.csv("../results/tailsDiabetesComplexmTAL.csv")
rownames(tails) <- NULL

colnames(iterprod) <- c("nums", "CI", "CIII", "CIV", "ATP", "Leak", "PO2")

diabetesTAL <- dplyr::inner_join(iterprod, tails, by = "nums")
diabetesTAL$V1 <- diabetesTAL$t <- diabetesTAL$X <- NULL
#setwd("..")
#setwd("..)
setwd("../liverScripts/modelScripts")

tails <- read.csv("../results/tailsDiabetes.csv")
rownames(tails) <- NULL

iterprod <- read.csv("../results/iterprodDiabetes.csv")

diabetesLiv <- dplyr::inner_join(iterprod, tails, by = "nums")
diabetesLiv$X.x <- diabetesLiv$X.y <- diabetesLiv$V1 <- diabetesLiv$t <- NULL

diabetesPT$PO2 <- as.factor(diabetesPT$PO2)
diabetesTAL$PO2 <- as.factor(diabetesTAL$PO2)
diabetesLiv$PO2 <- as.factor(diabetesLiv$PO2)

library(ggplot2)

g <- ggplot(diabetesPT, aes(x = 1000*ATP_c)) +
  geom_histogram(aes(fill = PO2)) +
  xlab("PT ATP (mM)") +
  ylab("Frequency") + 
  xlim(c(0, 3)) +
  geom_vline(xintercept = diabetesPT[348, ]$ATP_c*1000, col = "red") +
  theme(axis.title = element_text(size = 18), legend.position = "none", 
        axis.text = element_text(size = 18))

h <- ggplot(diabetesTAL, aes(x = 1000*ATP_c)) +
  geom_histogram(aes(fill = PO2)) +
  xlab("mTAL ATP (mM)") +
  xlim(c(0, 3)) +
  ylab("Frequency") + 
  geom_vline(xintercept = diabetesTAL[348, ]$ATP_c*1000, col = "red") +
  theme(axis.title = element_text(size = 18), legend.position = "none", 
        axis.text = element_text(size = 18))

i <- ggplot(diabetesLiv, aes(x = 1000*ATP_c)) +
  geom_histogram(aes(fill = PO2)) +
  xlab("Hepatocyte ATP (mM)") +
  xlim(c(0, 7.5)) +
  ylab("Frequency") + 
  geom_vline(xintercept = diabetesLiv[348, ]$ATP_c*1000, col = "red") +
  theme(axis.title = element_text(size = 18), legend.title = element_text(size = 18), 
        axis.text = element_text(size = 18), legend.text = element_text(size = 18))
i$labels$fill <- "Relative\nOxygen\nTension"

a <- ggplot(diabetesPT, aes(x = QH2_x/(0.00178*0.83))) +
  geom_histogram(aes(fill = PO2)) +
  xlab("PT Coenzyme Q Redox State") +
  xlim(c(0, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = diabetesPT[348, ]$QH2_x/(0.00178*0.83), col = "red") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18),
        legend.position = "none")

b <- ggplot(diabetesTAL, aes(x = QH2_x/(0.00178*0.83))) +
  geom_histogram(aes(fill = PO2)) +
  xlab("mTAL Coenzyme Q Redox State") +
  xlim(c(0, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = diabetesTAL[348, ]$QH2_x/(0.83*0.00178), col = "red") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18),
        legend.position = "none")

cg <- ggplot(diabetesLiv, aes(x = QH2_x/(6.49e-3*0.2))) +
  geom_histogram(aes(fill = PO2)) +
  xlab("Hepatocyte Coenzyme Q Redox State") +
  xlim(c(-0.025, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = diabetesLiv[423, ]$QH2_x/(6.49e-3*0.2), col = "red") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))
cg$labels$fill <- "Relative\nOxygen\nTension"

pdf("../dataVis/multipanelCoQDiabetes.pdf", width = 12)
gridExtra::grid.arrange(a, b, cg, nrow = 1, widths = c(1, 1, 1.2))
dev.off()

d <- ggplot(diabetesPT, aes(x = Cred_i/2.148e-3)) +
  geom_histogram(aes(fill = PO2)) +
  xlab("PT Cytochrome C Redox State") +
  xlim(c(0, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = diabetesPT[348, ]$Cred_i/2.148e-3, col = "red") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18),
        legend.position = "none")

e <- ggplot(diabetesTAL, aes(x = Cred_i/2.148e-3)) +
  geom_histogram(aes(fill = PO2)) +
  xlab("mTAL Cytochrome C Redox State") +
  xlim(c(0, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = diabetesTAL[348, ]$Cred_i/2.148e-3, col = "red") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18),
        legend.position = "none")

f <- ggplot(diabetesLiv, aes(x = Cred_i/0.0058387)) +
  geom_histogram(aes(fill = PO2)) +
  xlab("Hepatocyte Cytochrome C Redox State") +
  xlim(c(0, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = diabetesLiv[423, ]$Cred_i/0.0058387, col = "red") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))
f$labels$fill <- "Relative\nOxygen\nTension"

pdf("../dataVis/multipanelCytCDiabetes.pdf", width = 12)
gridExtra::grid.arrange(d, e, f, nrow = 1, widths = c(1, 1, 1.2))
dev.off()

## Figure 5.5
pdf("../dataVis/multipanelDiabetes.pdf", width = 14)
gridExtra::grid.arrange(g, h, i, a, b, cg, d, e, f, nrow = 3, widths = c(1, 1, 1.2))
dev.off()

## Multipanel with just uncoupling and hypoxia

olddiabetesPT <- diabetesPT
olddiabetesTAL <- diabetesTAL
olddiabetesLiv <- diabetesLiv

diabetesPT <- diabetesPT[diabetesPT$CI == 1.0 & diabetesPT$CIII == 1.0 & 
                           diabetesPT$CIV == 1.0 & diabetesPT$ATP == 1.0, ]
diabetesTAL <- diabetesTAL[diabetesTAL$CI == 1.0 & diabetesTAL$CIII == 1.0 & 
                             diabetesTAL$CIV == 1.0 & diabetesTAL$ATP == 1.0, ]
diabetesLiv <- diabetesLiv[diabetesLiv$CI == 1.0 & diabetesLiv$CIII == 1.0 & 
                             diabetesLiv$CIV == 1.0 & diabetesLiv$ATP == 1.0, ]

g <- ggplot(diabetesPT, aes(x = 1000*ATP_c)) +
  geom_histogram(aes(fill = PO2)) +
  xlab("PT ATP (mM)") +
  ylab("Frequency") + 
  xlim(c(0, 3)) +
  geom_vline(xintercept = diabetesPT[3, ]$ATP_c*1000, col = "red") +
  theme(axis.title = element_text(size = 18), legend.position = "none", 
        axis.text = element_text(size = 18))

h <- ggplot(diabetesTAL, aes(x = 1000*ATP_c)) +
  geom_histogram(aes(fill = PO2)) +
  xlab("mTAL ATP (mM)") +
  xlim(c(0, 3)) +
  ylab("Frequency") + 
  geom_vline(xintercept = diabetesTAL[3, ]$ATP_c*1000, col = "red") +
  theme(axis.title = element_text(size = 18), legend.position = "none", 
        axis.text = element_text(size = 18))

i <- ggplot(diabetesLiv, aes(x = 1000*ATP_c)) +
  geom_histogram(aes(fill = PO2)) +
  xlab("Hepatocyte ATP (mM)") +
  xlim(c(0, 8)) +
  ylab("Frequency") + 
  geom_vline(xintercept = diabetesLiv[3, ]$ATP_c*1000, col = "red") +
  theme(axis.title = element_text(size = 18), legend.title = element_text(size = 18), 
        axis.text = element_text(size = 18), legend.text = element_text(size = 18))
i$labels$fill <- "Relative\nOxygen\nTension"

pdf("../dataVis/multipanelUncouplingHypoxia.pdf")
gridExtra::grid.arrange(g, h, i)
dev.off()

a <- ggplot(diabetesPT, aes(x = QH2_x/(0.00178*0.83))) +
  geom_histogram(aes(fill = PO2)) +
  xlab("PT Coenzyme Q Redox State") +
  xlim(c(0, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = diabetesPT[3, ]$QH2_x/(0.00178*0.83), col = "red") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18),
        legend.position = "none")

b <- ggplot(diabetesTAL, aes(x = QH2_x/(0.00178*0.83))) +
  geom_histogram(aes(fill = PO2)) +
  xlab("mTAL Coenzyme Q Redox State") +
  xlim(c(0, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = diabetesTAL[3, ]$QH2_x/(0.83*0.00178), col = "red") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18),
        legend.position = "none")

cg <- ggplot(diabetesLiv, aes(x = QH2_x/(6.49e-3*0.2))) +
  geom_histogram(aes(fill = PO2)) +
  xlab("Hepatocyte Coenzyme Q Redox State") +
  xlim(c(-0.025, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = diabetesLiv[3, ]$QH2_x/(6.49e-3*0.2), col = "red") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))
cg$labels$fill <- "Relative\nOxygen\nTension"

d <- ggplot(diabetesPT, aes(x = Cred_i/2.148e-3)) +
  geom_histogram(aes(fill = PO2)) +
  xlab("PT Cytochrome C Redox State") +
  xlim(c(0, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = diabetesPT[3, ]$Cred_i/2.148e-3, col = "red") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18),
        legend.position = "none")

e <- ggplot(diabetesTAL, aes(x = Cred_i/2.148e-3)) +
  geom_histogram(aes(fill = PO2)) +
  xlab("mTAL Cytochrome C Redox State") +
  xlim(c(0, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = diabetesTAL[3, ]$Cred_i/2.148e-3, col = "red") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18),
        legend.position = "none")

f <- ggplot(diabetesLiv, aes(x = Cred_i/0.0058387)) +
  geom_histogram(aes(fill = PO2)) +
  xlab("Hepatocyte Cytochrome C Redox State") +
  xlim(c(0, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = diabetesLiv[3, ]$Cred_i/0.0058387, col = "red") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))
f$labels$fill <- "Relative\nOxygen\nTension"

## Figure 5.4
pdf("../dataVis/multipanelUncouplingHypoxia2.pdf", width = 14)
gridExtra::grid.arrange(g, h, i, a, b, cg, d, e, f, nrow = 3, widths = c(1, 1, 1.2))
dev.off()

diabetesPT <- olddiabetesPT[olddiabetesPT$Leak == 1. & olddiabetesPT$PO2 == 1., ]
diabetesTAL <- olddiabetesTAL[olddiabetesTAL$Leak == 1. & olddiabetesTAL$PO2 == 1., ]
diabetesLiv <- olddiabetesLiv[olddiabetesLiv$Leak == 1. & olddiabetesLiv$PO2 == 1., ]

diabetesPT$CI <- as.factor(diabetesPT$CI)
diabetesTAL$CI <- as.factor(diabetesTAL$CI)
diabetesLiv$CI <- as.factor(diabetesLiv$CI)

diabetesPT$CIII <- as.factor(diabetesPT$CIII)
diabetesTAL$CIII <- as.factor(diabetesTAL$CIII)
diabetesLiv$CIII <- as.factor(diabetesLiv$CIII)

diabetesPT$CIV <- as.factor(diabetesPT$CIV)
diabetesTAL$CIV <- as.factor(diabetesTAL$CIV)
diabetesLiv$CIV <- as.factor(diabetesLiv$CIV)

a <- ggplot(diabetesPT, aes(x = QH2_x/(0.00178*0.83))) +
  geom_histogram(aes(fill = CI)) +
  xlab("PT Coenzyme Q Redox State") +
  xlim(c(0, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = diabetesPT[348, ]$QH2_x/(0.00178*0.83), col = "red") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))
a$labels$fill <- "Complex I\nActivity"

b <- ggplot(diabetesTAL, aes(x = QH2_x/(0.00178*0.83))) +
  geom_histogram(aes(fill = CI)) +
  xlab("mTAL Coenzyme Q Redox State") +
  xlim(c(0, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = diabetesTAL[348, ]$QH2_x/(0.83*0.00178), col = "red") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))
b$labels$fill <- "Complex I\nActivity"

cg <- ggplot(diabetesLiv, aes(x = QH2_x/(6.49e-3*0.2))) +
  geom_histogram(aes(fill = CI)) +
  xlab("Hepatocyte Coenzyme Q Redox State") +
  xlim(c(-0.025, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = diabetesLiv[423, ]$QH2_x/(6.49e-3*0.2), col = "red") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) +
  scale_fill_manual(values = gg_color_hue(4)[c(2, 3, 4)])
cg$labels$fill <- "Complex I\nActivity"

## Figure 5.2
pdf("../dataVis/multipanelQH2OXPHOS.pdf")
gridExtra::grid.arrange(a, b, cg)
dev.off()

g <- ggplot(diabetesPT, aes(x = 1000*ATP_c)) +
  geom_histogram(aes(fill = CIV)) +
  xlab("PT ATP (mM)") +
  ylab("Frequency") + 
  xlim(c(0, 3)) +
  theme(axis.title = element_text(size = 18), legend.title = element_text(size = 18), 
        axis.text = element_text(size = 18), legend.text = element_text(size = 18))
g$labels$fill <- "Complex IV\nActivity"

h <- ggplot(diabetesTAL, aes(x = 1000*ATP_c)) +
  geom_histogram(aes(fill = CIV)) +
  xlab("mTAL ATP (mM)") +
  xlim(c(0, 3)) +
  ylab("Frequency") + 
  theme(axis.title = element_text(size = 18), legend.title = element_text(size = 18), 
        axis.text = element_text(size = 18), legend.text = element_text(size = 18))
h$labels$fill <- "Complex IV\nActivity"

i <- ggplot(diabetesLiv, aes(x = 1000*ATP_c)) +
  geom_histogram(aes(fill = CIV)) +
  xlab("Hepatocyte ATP (mM)") +
  xlim(c(0, 7.5)) +
  ylab("Frequency") + 
  theme(axis.title = element_text(size = 18), legend.title = element_text(size = 18), 
        axis.text = element_text(size = 18), legend.text = element_text(size = 18)) +
  scale_fill_manual(values = c(gg_color_hue(4), gg_color_hue(5)[5]))
i$labels$fill <- "Complex IV\nActivity"

## Figure 5.1
pdf("../dataVis/multipanelATPoxphosDiabetes.pdf")
gridExtra::grid.arrange(g, h, i)
dev.off()


d <- ggplot(diabetesPT, aes(x = Cred_i/2.148e-3)) +
  geom_histogram(aes(fill = CIV)) +
  xlab("PT Cytochrome C Redox State") +
  xlim(c(0, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = diabetesPT[348, ]$Cred_i/2.148e-3, col = "red") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))
d$labels$fill <- "Complex IV\nActivity"

e <- ggplot(diabetesTAL, aes(x = Cred_i/2.148e-3)) +
  geom_histogram(aes(fill = CIV)) +
  xlab("mTAL Cytochrome C Redox State") +
  xlim(c(0, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = diabetesTAL[348, ]$Cred_i/2.148e-3, col = "red") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))
e$labels$fill <- "Complex IV\nActivity"

f <- ggplot(diabetesLiv, aes(x = Cred_i/0.0058387)) +
  geom_histogram(aes(fill = CIV)) +
  xlab("Hepatocyte Cytochrome C Redox State") +
  xlim(c(0, 1)) +
  ylab("Frequency") +
  geom_vline(xintercept = diabetesLiv[423, ]$Cred_i/0.0058387, col = "red") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  scale_fill_manual(values = c(gg_color_hue(4), gg_color_hue(5)[5]))
f$labels$fill <- "Complex IV\nActivity"

## Figure 5.3
pdf("../dataVis/multipanelCytCOXPHOSdiabetes.pdf")
gridExtra::grid.arrange(d, e, f)
dev.off()

