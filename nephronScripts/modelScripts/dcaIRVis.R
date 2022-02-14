if (!grepl("modelScripts", getwd())) {
        setwd("./modelScripts")
}
reperfusionDCA <- read.csv("../results/resultsReperfusionPoolDCA.csv")
reperfusionPYR <- read.csv("../results/resultsReperfusionPoolPYR.csv")
reperfusionNoDCA <- read.csv("../results/resultsReperfusionPool12.csv")

pdf('../dataVis/dcaIR.pdf', width = 15)
par(mfrow = c(1, 2), cex.axis = 2, cex.lab = 2.5, mar = c(5.1, 5, 4.1, 3))
plot(shortreperfusionDCA$t, shortreperfusionDCA$NADH_x/0.824e-3, cex = 0.1, col = "blue",
     ylim = c(0, 1), ylab = "Redox State of NADH", xlab = "Time (s)")
lines(shortreperfusionDCA$t, shortreperfusionDCA$NADH_x/0.824e-3, cex = 0.1, col = "blue")
lines(shortreperfusionNoDCA$t, shortreperfusionNoDCA$NADH_x/0.824e-3, col = "red")

# plot(shortreperfusionDCA$t, shortreperfusionDCA$QH2_x/6.49e-3, cex = 0.1, col = "blue",
#      ylim = c(0, 1), ylab = "Redox State of CoQ", xlab = "Time (s)")
# lines(shortreperfusionDCA$t, shortreperfusionDCA$QH2_x/6.49e-3, cex = 0.1, col = "blue")
# lines(shortreperfusionNoDCA$t, shortreperfusionNoDCA$QH2_x/6.49e-3, cex = 0.1, col = "red")

dpH <- -0.059*log(shortreperfusionDCA$H_x/shortreperfusionDCA$H_i)*1000
dpHNoDCA <-  -0.059*log(shortreperfusionNoDCA$H_x/shortreperfusionNoDCA$H_i)*1000
plot(shortreperfusionDCA$t, shortreperfusionDCA$dPsi+dpH, cex = 0.1, col = "blue", ylim = c(150, 200),
     ylab = "Proton Motive Force (mV)", xlab = "Time (s)")
lines(shortreperfusionDCA$t, shortreperfusionDCA$dPsi+dpH, cex = 0.1, col = "blue")
lines(shortreperfusionNoDCA$t, shortreperfusionNoDCA$dPsi+dpHNoDCA, col = "red")

# plot(shortreperfusionDCA$t, shortreperfusionDCA$SUC_x*1000, cex = 0.1, col = "blue", ylim = c(0, 40),
#      ylab = "Succinate (mM)", xlab = "Time (s)")
# lines(shortreperfusionDCA$t, shortreperfusionDCA$SUC_x*1000, cex = 0.1, col = "blue")
# lines(shortreperfusionNoDCA$t, shortreperfusionNoDCA$SUC_x*1000, col = "red")
par(mfrow = c(1, 1))
dev.off()

## Pyruvate case
shortreperfusionDCA <- reperfusionPYR[reperfusionPYR$t < 250, ]

pdf('../dataVis/pyrIR.pdf', width = 15)
par(mfrow = c(1, 4), cex.axis = 2, cex.lab = 2.5, mar = c(5.1, 5, 4.1, 3))
plot(shortreperfusionDCA$t, shortreperfusionDCA$NADH_x/4.85e-3, cex = 0.1, col = "blue",
     ylim = c(0, 1), ylab = "Redox State of NADH", xlab = "Time (s)")
lines(shortreperfusionDCA$t, shortreperfusionDCA$NADH_x/4.85e-3, cex = 0.1, col = "blue")
lines(shortreperfusionNoDCA$t, shortreperfusionNoDCA$NADH_x/4.85e-3, col = "red")

plot(shortreperfusionDCA$t, shortreperfusionDCA$QH2_x/6.49e-3, cex = 0.1, col = "blue",
     ylim = c(0, 1), ylab = "Redox State of CoQ", xlab = "Time (s)")
lines(shortreperfusionDCA$t, shortreperfusionDCA$QH2_x/6.49e-3, cex = 0.1, col = "blue")
lines(shortreperfusionNoDCA$t, shortreperfusionNoDCA$QH2_x/6.49e-3, cex = 0.1, col = "red")

dpH <- -0.059*log(shortreperfusionDCA$H_x/shortreperfusionDCA$H_i)*1000
dpHNoDCA <-  -0.059*log(shortreperfusionNoDCA$H_x/shortreperfusionNoDCA$H_i)*1000
plot(shortreperfusionDCA$t, shortreperfusionDCA$dPsi+dpH, cex = 0.1, col = "blue", ylim = c(150, 230),
     ylab = "Proton Motive Force (mV)", xlab = "Time (s)")
lines(shortreperfusionDCA$t, shortreperfusionDCA$dPsi+dpH, cex = 0.1, col = "blue")
lines(shortreperfusionNoDCA$t, shortreperfusionNoDCA$dPsi+dpHNoDCA, col = "red")

plot(shortreperfusionDCA$t, shortreperfusionDCA$SUC_x*1000, cex = 0.1, col = "blue", ylim = c(0, 40),
     ylab = "Succinate (mM)", xlab = "Time (s)")
lines(shortreperfusionDCA$t, shortreperfusionDCA$SUC_x*1000, cex = 0.1, col = "blue")
lines(shortreperfusionNoDCA$t, shortreperfusionNoDCA$SUC_x*1000, col = "red")
par(mfrow = c(1, 1))
dev.off()