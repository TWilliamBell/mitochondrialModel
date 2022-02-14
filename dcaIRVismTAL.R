if (!grepl("modelScripts", getwd())) {
  setwd("modelScripts")
}

reperfusionDCA <- data.table::fread("../results/resultsReperfusionPoolDCAmTAL.csv")
reperfusionPYR <- data.table::fread("../results/resultsReperfusionPoolPYRmTAL.csv")
reperfusionNoDCA <- feather::read_feather("../results/resultsReperfusionPool12mTAL.feather")

par(mfrow = c(1, 1))
plot(reperfusionDCA$t, reperfusionDCA$NADH_x/0.824e-3, cex = 0.1, col = "red",
     ylim = c(0, 1), ylab = "Redox State of NADH", xlab = "Time (s)")
lines(reperfusionDCA$t, reperfusionDCA$NADH_x/0.824e-3, cex = 0.1, col = "red")

plot(reperfusionDCA$t, reperfusionDCA$QH2_x/2.148e-3, cex = 0.1, col = "red",
     ylim = c(0, 1))
lines(reperfusionDCA$t, reperfusionDCA$QH2_x/2.148e-3, cex = 0.1, col = "red")

plot(reperfusionDCA$t, reperfusionDCA$Cred_i, cex = 0.1, col = "red")
lines(reperfusionDCA$t, reperfusionDCA$Cred_i, cex = 0.1, col = "red")

plot(reperfusionDCA$t, reperfusionDCA$dPsi, cex = 0.1, col = "red")
lines(reperfusionDCA$t, reperfusionDCA$dPsi, cex = 0.1, col = "red")

shortreperfusionDCA <- reperfusionDCA[reperfusionDCA$t < 50, ]
shortreperfusionNoDCA <- reperfusionNoDCA[reperfusionNoDCA$t < 50, ]

pdf('../dataVis/dcaIRmTAL.pdf', width = 15)
par(mfrow = c(1, 3), cex.axis = 2, cex.lab = 2.5, mar = c(5.1, 5, 4.1, 3))
plot(shortreperfusionDCA$t, shortreperfusionDCA$NADH_x/0.824e-3, cex = 0.1, col = "blue",
     ylim = c(0, 1), ylab = "Redox State of NADH", xlab = "Time (s)")
lines(shortreperfusionDCA$t, shortreperfusionDCA$NADH_x/0.824e-3, cex = 0.1, col = "blue")
lines(shortreperfusionNoDCA$t, shortreperfusionNoDCA$NADH_x/0.824e-3, col = "red")

plot(shortreperfusionDCA$t, shortreperfusionDCA$QH2_x/2.148e-3, cex = 0.1, col = "blue",
     ylim = c(0, 1), ylab = "Redox State of CoQ", xlab = "Time (s)")
lines(shortreperfusionDCA$t, shortreperfusionDCA$QH2_x/2.148e-3, cex = 0.1, col = "blue")
lines(shortreperfusionNoDCA$t, shortreperfusionNoDCA$QH2_x/2.148e-3, cex = 0.1, col = "red")

dpH <- -0.059*log(shortreperfusionDCA$H_x/shortreperfusionDCA$H_i)*1000
dpHNoDCA <-  -0.059*log(shortreperfusionNoDCA$H_x/shortreperfusionNoDCA$H_i)*1000
plot(shortreperfusionDCA$t, shortreperfusionDCA$dPsi+dpH, cex = 0.1, col = "blue", ylim = c(150, 200),
     ylab = "Proton Motive Force (mV)", xlab = "Time (s)")
lines(shortreperfusionDCA$t, shortreperfusionDCA$dPsi+dpH, cex = 0.1, col = "blue")
lines(shortreperfusionNoDCA$t, shortreperfusionNoDCA$dPsi+dpHNoDCA, col = "red")

par(mfrow = c(1, 1))
dev.off()

## Pyruvate case
shortreperfusionDCA <- reperfusionPYR[reperfusionPYR$t < 50, ]

pdf('../dataVis/pyrIRmTAL.pdf', width = 15)
par(mfrow = c(1, 3), cex.axis = 2, cex.lab = 2.5, mar = c(5.1, 5, 4.1, 3))
plot(shortreperfusionDCA$t, shortreperfusionDCA$NADH_x/0.824e-4, cex = 0.1, col = "blue",
     ylim = c(0, 1), ylab = "Redox State of NADH", xlab = "Time (s)")
lines(shortreperfusionDCA$t, shortreperfusionDCA$NADH_x/0.824e-4, cex = 0.1, col = "blue")
lines(shortreperfusionNoDCA$t, shortreperfusionNoDCA$NADH_x/0.824e-4, col = "red")

plot(shortreperfusionDCA$t, shortreperfusionDCA$QH2_x/2.148e-3, cex = 0.1, col = "blue",
     ylim = c(0, 1), ylab = "Redox State of CoQ", xlab = "Time (s)")
lines(shortreperfusionDCA$t, shortreperfusionDCA$QH2_x/2.148e-3, cex = 0.1, col = "blue")
lines(shortreperfusionNoDCA$t, shortreperfusionNoDCA$QH2_x/2.148e-3, cex = 0.1, col = "red")

dpH <- -0.059*log(shortreperfusionDCA$H_x/shortreperfusionDCA$H_i)*1000
dpHNoDCA <-  -0.059*log(shortreperfusionNoDCA$H_x/shortreperfusionNoDCA$H_i)*1000
plot(shortreperfusionDCA$t, shortreperfusionDCA$dPsi+dpH, cex = 0.1, col = "blue", ylim = c(150, 200),
     ylab = "Proton Motive Force (mV)", xlab = "Time (s)")
lines(shortreperfusionDCA$t, shortreperfusionDCA$dPsi+dpH, cex = 0.1, col = "blue")
lines(shortreperfusionNoDCA$t, shortreperfusionNoDCA$dPsi+dpHNoDCA, col = "red")

par(mfrow = c(1, 1))
dev.off()