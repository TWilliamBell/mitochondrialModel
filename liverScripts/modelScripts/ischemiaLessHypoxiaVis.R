if (!grepl("modelScripts", getwd())) {
  setwd("liverScripts/modelScripts")
}

reperfusion <- data.table::fread("../results/resultsReperfusionPoolLessHypox.csv")
reperfusion <- rbind(reperfusion[1, ], reperfusion)
reperfusion$t[1] <- -50
shortreperfusion <- reperfusion[reperfusion$t < 200, ]

pdf("../dataVis/ischemiaLessHypoxia.pdf", width = 15)
par(mfrow = c(1, 4), cex.axis = 2, cex.lab = 2.5, mar = c(5.1, 5, 4.1, 3))

dpH <- -0.059*log(shortreperfusion$H_x/shortreperfusion$H_i)*1000

plot(shortreperfusion$t, shortreperfusion$NADH_x/4.85e-3, ylim = c(0, 1), cex = 0.1, col = "red",
     xlab = "Time (s)", ylab = "Proportion of the NADH/NAD+ Pool in Reduced State")
lines(shortreperfusion$t, shortreperfusion$NADH_x/4.85e-3, col = "red")

plot(shortreperfusion$t, shortreperfusion$QH2_x/6.49e-3, ylim = c(0, 1), cex = 0.1, col = "red",
     xlab = "Time (s)", ylab = "Proportion of the Coenzyme Q Pool in Reduced State")

plot(shortreperfusion$t, shortreperfusion$dPsi+dpH, cex = 0.1, col = "red", xlab = "Time (s)", 
     ylab = "Proton Motive Force (mV)")

plot(shortreperfusion$t, 1000*shortreperfusion$SUC_x, cex = 0.1, col = "red", xlab = "Time (s)",
     ylab = "Matrix Succinate (mM)")
lines(shortreperfusion$t, 1000*shortreperfusion$SUC_x, col = "red")
dev.off()

## 
reperfusionPoolOXPHOS12 <- data.table::fread(
  "../results/resultsReperfusionOXPHOSPool12.csv")
reperfusionPoolOXPHOS12 <- rbind(reperfusionPoolOXPHOS12[1, ], reperfusionPoolOXPHOS12)
reperfusionPoolOXPHOS12$t[1] <- -50
reperfusionPool12 <- data.table::fread(
  "../results/resultsReperfusionPool12.csv")
reperfusionPool12 <- rbind(reperfusionPool12[1, ], reperfusionPool12)
reperfusionPool12$t[1] <- -50

shortReperfusionPoolOXPHOS12 <- reperfusionPoolOXPHOS12[
reperfusionPoolOXPHOS12$t < 200, ]
dpHPOX <- -0.059*log(shortReperfusionPoolOXPHOS12$H_x/shortReperfusionPoolOXPHOS12$H_i)*1000

shortReperfusionPool12 <- reperfusionPool12[
  reperfusionPool12$t < 200, ]
dpHP <- -0.059*log(shortReperfusionPool12$H_x/shortReperfusionPool12$H_i)*1000

pdf("../dataVis/bigfnIschemia.pdf", width = 15)
par(mfrow = c(1, 4), cex.axis = 2, cex.lab = 2.5, mar = c(5.1, 5, 4.1, 3))
plot(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$NADH_x/4.85e-3,
     col = "orange", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
     ylab = "Proportion of NADH/NAD+ Pool in Reduced State")
lines(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$NADH_x/4.85e-3,
      col = "orange")
lines(shortReperfusionPool12$t, shortReperfusionPool12$NADH_x/4.85e-3,
      col = "red")
lines(shortreperfusion$t, shortreperfusion$NADH_x/4.85e-3, col = "blue")
plot(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$QH2_x/6.49e-3,
     col = "orange", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
     ylab = "Proportion of Coenzyme Q Pool in Reduced State")
lines(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$QH2_x/6.49e-3,
      col = "orange")
lines(shortreperfusion$t, shortreperfusion$QH2_x/6.49e-3,
      col = "blue")
lines(shortReperfusionPool12$t, shortReperfusionPool12$QH2_x/6.49e-3,
      col = "red")
plot(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$dPsi+dpHPOX,
     col = "orange", cex = 0.1, xlab = "Time (s)", ylim = c(100, 250),
     ylab = "Proton Motive Force (mV)")
lines(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$dPsi+dpHPOX,
      col = "orange")
lines(shortreperfusion$t, shortreperfusion$dPsi+dpH,
      col = "blue")
lines(shortReperfusionPool12$t, shortReperfusionPool12$dPsi+dpHP,
      col = "red")
plot(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$SUC_x*1000, cex = 0.1,
     col = "orange", xlab = "Time (s)", ylab = "Matrix Succinate (mM)")
lines(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$SUC_x*1000,
      col = "orange")
lines(shortreperfusion$t, 1000*shortreperfusion$SUC_x, col = "blue")
lines(shortReperfusionPool12$t, shortReperfusionPool12$SUC_x*1000,
      col = "red")
dev.off()
