## In order to run this file, you should first run fuelRich.py, hypoxia.py, 
## and mainDynamical.py

if (!grepl("liverScripts/modelScripts", getwd())) {
  setwd("./liverScripts/modelScripts")
}

default <- data.table::fread("../results/resultsATP.csv")
fuel1 <- data.table::fread("../results/resultsFuel0.csv")
fuel2 <- data.table::fread("../results/resultsFuel1.csv")

hypoxiaExtreme <- data.table::fread("../results/resultsHypoxiaExtreme.csv")

plot(hypoxiaExtreme$t[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000], 
     hypoxiaExtreme$ATP_c[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000],
     cex = 0.1, col = 'red', xlab = "Time (s)", ylab = "[ATP]_c (M)", ylim = c(0, 0.007))
lines(hypoxiaExtreme$t[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000], 
      hypoxiaExtreme$ATP_c[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000], col = "red")


plot(hypoxiaExtreme$t[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000], 
     hypoxiaExtreme$Cred_i[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000]/4.39e-3,
     cex = 0.1, col = 'red', xlab = "Time (s)", ylab = "Cytochrome Reduction State", 
     ylim = c(0, 1))
lines(hypoxiaExtreme$t[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000], 
      hypoxiaExtreme$Cred_i[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000]/4.39e-3, 
      col = "red")

# plot(hypoxiaExtreme$t[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000],
#      hypoxiaExtreme$NADH_x[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000]/9.4e-3,
#      cex = 0.1, col = 'red', xlab = "Time (s)", ylab = "NADH Reduction State",
#      ylim = c(0, 1))
# lines(hypoxiaExtreme$t[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000], 
#       hypoxiaExtreme$NADH_x[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000]/9.4e-3, 
#       col = "red")
# 
# plot(hypoxiaExtreme$t[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000],
#      hypoxiaExtreme$QH2_x[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000]/6.49e-3,
#      cex = 0.1, col = 'red', xlab = "Time (s)", ylab = "CoQ Reduction State",
#      ylim = c(0, 1))
# lines(hypoxiaExtreme$t[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000], 
#       hypoxiaExtreme$QH2_x[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000]/6.49e-3, 
#       col = "red")

# plot(hypoxiaExtreme$t[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000], 
#      hypoxiaExtreme$dPsi[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000],
#      cex = 0.1, col = 'red', xlab = "Time (s)", ylab = "Electrical Potential Gradient")
# lines(hypoxiaExtreme$t[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000], 
#       hypoxiaExtreme$dPsi[hypoxiaExtreme$t < 25000 & hypoxiaExtreme$t > 5000], 
#       col = "red")

tailDefault <- tail(default, 1)
tailDefault$V1 <- tailDefault$t <- NULL
tailFuel1 <- tail(fuel1, 1)
tailFuel1$V1 <- tailFuel1$t <- NULL
tailFuel2 <- tail(fuel2, 1)
tailFuel2$V1 <- tailFuel2$t <- NULL

CRed <- c(tailDefault$Cred_i, tailFuel1$Cred_i)/4.39e-3  #, tailFuel2$Cred_i)
QRed <- c(tailDefault$QH2_x, tailFuel1$QH2_x)/6.49e-3  #, tailFuel2$QH2_x)
NADH <- c(tailDefault$NADH_x, tailFuel1$NADH_x)/9.4e-3  #, tailFuel2$NADH_x)
dPsi <- c(tailDefault$dPsi, tailFuel1$dPsi)
pmf <- c(tailDefault$dPsi, tailFuel1$dPsi)+c(-0.059*log(tailDefault$H_x/tailDefault$H_i)*1000,
                                             -0.059*log(tailFuel1$H_x/tailFuel1$H_i)*1000)

## Figure 4.3
pdf("../dataVis/resultsGlycolysis.pdf")
par(cex.axis = 1.5, cex.lab = 1.5)
#par(mfrow = c(1, 2), mar = c(5, 4, 4, 0)+0.1, cex.axis = 2, cex.lab = 1.4)
#barplot(CRed, names.arg = c("1x", "4x"), col = "mistyrose2", xlab = "CytC")
#barplot(QRed, names.arg = c("1x", "4x"), col = "mistyrose2", xlab = "CoQ")
#barplot(NADH, names.arg = c("1x", "4x"), col = "mistyrose2", xlab = "NADH Reduction State",
#        ylim = c(0, 1))
#barplot(dPsi, names.arg = c("1x", "4x"), col = "mistyrose2", xlab = "Potential Gradient",
#        xpd = F, ylim = c(150, 210))
barplot(pmf, names.arg = c("1x", "4x"), col = "mistyrose2", ylab = "Proton Motive Force (mV)",
        xlab = "Fold-Change in Glycolytic Activity",
        xpd = F, ylim = c(0, 250))
#par(mfrow = c(1, 1))
dev.off()

CRed <- c(tailDefault$Cred_i, tailFuel2$Cred_i)
QRed <- c(tailDefault$QH2_x, tailFuel2$QH2_x)
NADH <- c(tailDefault$NADH_x, tailFuel2$NADH_x)
dPsi <- c(tailDefault$dPsi, tailFuel2$dPsi)

pdf("../dataVis/resultsHypoxicROS.pdf", width = 12)
par(mfrow = c(1, 2), mar = c(5, 6, 4, 0)+0.1, cex.axis = 2, cex.lab = 2)
#barplot(CRed, names.arg = c("36 mmHg", "5 mmHg"), col = "mistyrose2")#, ylim = c(0, 0.002))
barplot(NADH, names.arg = c("36 mmHg", "5 mmHg"), col = "mistyrose2", ylim = c(0, 1),
        ylab = "NADH Reduction State")
barplot(dPsi, names.arg = c("36 mmHg", "5 mmHg"), col = "mistyrose2", ylim = c(150, 200),
        xpd = F, ylab = "Electrical Potential Gradient")
par(mfrow = c(1, 1))
dev.off()
