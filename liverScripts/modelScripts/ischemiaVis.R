## In order to run this file, you should first run ischemia.py

if (!grepl("liverScripts/modelScripts", getwd())) {
  setwd("./liverScripts/modelScripts")
}

ischemia <- data.table::fread("../results/resultsIschemia.csv")

baseline <- data.table::fread("../results/resultsATP.csv")
# reperfusion1 <- data.table::fread("../results/resultsReperfusion1.csv")
# reperfusion2 <- data.table::fread("../results/resultsReperfusion2.csv")
# reperfusion12 <- data.table::fread("../results/resultsReperfusion12.csv")

#reperfusionOXPHOS <- feather::read_feather("../results/resultsReperfusionOXPHOS1.feather")
#reperfusionOXPHOS2 <- feather::read_feather("../results/resultsReperfusionOXPHOS2.feather")
#reperfusionOXPHOS12 <- feather::read_feather("../results/resultsReperfusionOXPHOS12.feather")

#reperfusionPool1 <- feather::read_feather(
#  "../results/resultsReperfusionPool1.feather")
reperfusionPool12 <- data.table::fread(
  "../results/resultsReperfusionPool12.csv")

# reperfusionPoolOXPHOS1 <- data.table::fread(
#   "../results/resultsReperfusionOXPHOSPool1.csv")
# reperfusionPoolOXPHOS2 <- data.table::fread(
#   "../results/resultsReperfusionOXPHOSPool2.csv")
reperfusionPoolOXPHOS12 <- data.table::fread(
  "../results/resultsReperfusionOXPHOSPool12.csv")

# shortreperfusion1 <- reperfusion1[reperfusion1$t < 250, ]
# shortreperfusion2 <- reperfusion2[reperfusion2$t < 250, ]
# shortreperfusion12 <- reperfusion12[reperfusion12$t < 250, ]
# 
# pdf("../dataVis/ischemiaInjury.pdf")
# par(mfrow = c(1, 2), mar = c(5, 6, 4, 1.5)+0.1, cex.axis = 2, cex.lab = 2.5)
# plot(shortreperfusion1$t, shortreperfusion1$NADH_x/9.4e-3, cex = 0.1, 
#      ylim = c(0, 1), xlab = "Time (s)", ylab = "NADH Reduction State", col = "red")
# par(mar = c(5, 6, 4, 0)+0.1)
# plot(shortreperfusion1$t, shortreperfusion1$dPsi, cex = 0.1, 
#      ylim = c(150, 200), xlab = "Time (s)", ylab = "Electrical Potential Gradient", 
#      col = "red")
# par(mfrow = c(1, 1), mar = c(5, 4, 4, 2)+0.1)
# dev.off()

# # par(mfrow = c(1, 3), mar = c(5, 4, 4, 0) + 0.1)
# plot(shortreperfusion1$t, (shortreperfusion1$Cred_i/4.39e-3)*5, cex = 0.1, 
#      ylim = c(0, 1), xlab = "Time (s)", ylab = "CytC Reduction State", col = "red")
# plot(shortreperfusion1$t, shortreperfusion1$QH2_x/6.49e-3, cex = 0.1, 
#      ylim = c(0, 1), xlab = "Time (s)", ylab = "CoQ Reduction State", col = "red")
# 
# # par(mar = c(5, 2, 4, 1.5)+0.1)
# # plot(shortreperfusion2$t, shortreperfusion2$NADH_x, cex = 0.1, 
# #      ylim = c(0, 0.003), xlab = "Time (s)", ylab = "", col = "red")
# # par(mar = c(5, 1, 4, 2.5)+0.1)
# # plot(shortreperfusion12$t, shortreperfusion12$NADH_x, cex = 0.1, 
# #      ylim = c(0, 0.003), xlab = "Time (s)", ylab = "", col = "red")
# # par(mfrow = c(1, 1), mar = c(5, 4, 4, 2)+0.1)



# par(mfrow = c(1, 3), mar = c(5, 4, 4, 0) + 0.1)
# plot(shortreperfusion1$t, shortreperfusion1$QH2_x, cex = 0.1, 
#      ylim = c(0, 0.0008), xlab = "Time (s)", ylab = "QH_2 (M)", col = "red")
# par(mar = c(5, 2, 4, 1.5)+0.1)
# plot(shortreperfusion2$t, shortreperfusion2$QH2_x, cex = 0.1, col = "red",
#      ylim = c(0, 0.0008), xlab = "Time (s)", ylab = "")
# par(mar = c(5, 1, 4, 2.5)+0.1)
# plot(shortreperfusion12$t, shortreperfusion12$QH2_x, cex = 0.1, col = "red",
#      ylim = c(0, 0.0008), xlab = "Time (s)", ylab = "")
# par(mfrow = c(1, 1), mar = c(5, 4, 4, 2)+0.1)

# par(mfrow = c(1, 3))
# plot(shortreperfusion1$t, shortreperfusion1$Cred_i, cex = 0.1, 
#      ylim = c(0, 0.0015), xlab = "Time (s)", ylab = "cytochrome C (M)")
# plot(shortreperfusion2$t, shortreperfusion2$Cred_i, cex = 0.1, 
#      ylim = c(0, 0.0015), xlab = "Time (s)", ylab = "")
# plot(shortreperfusion12$t, shortreperfusion12$Cred_i, cex = 0.1, 
#      ylim = c(0, 0.0015), xlab = "Time (s)", ylab = "")
# par(mfrow = c(1, 1))

# shortreperfusionOXPHOS <- reperfusionOXPHOS[reperfusionOXPHOS$t < 200, ]
# shortreperfusionOXPHOS2 <- reperfusionOXPHOS2[reperfusionOXPHOS2$t < 200, ]
# 
# dpH <- -0.059*log(shortreperfusionOXPHOS$H_x/shortreperfusionOXPHOS$H_i)*1000
# 
# par(mfrow = c(1, 2))
# plot(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$NADH_x/9.4e-3,
#      xlab = "Time (s)", ylab = "NADH Proportion of Nicotinamide adenine Pool",
#      cex = 0.1, ylim = c(0, 1), col = "red")
# lines(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$NADH_x/9.4e-3, col = "red")
# plot(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$QH2_x/6.49e-3,
#      xlab = "Time (s)", ylab = "Proportion of Coenzyme Q in a Reduced State",
#      cex = 0.1, ylim = c(0, 1), col = "red")
# lines(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$QH2_x/6.49e-3, col = "red")
# par(mfrow = c(1, 1))
# # plot(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$dPsi,
# #      xlab = "Time (s)", ylab = "Electrical Potential Gradient (mV)",
# #      cex = 0.1, col = "red", ylim = c(140, 180))
# # lines(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$dPsi, col = "red")
# # plot(shortreperfusionOXPHOS$t, dpH, cex = 0.1, xlab = "Time (s)", 
# #      ylab= "pH Difference (mV)", col = "red", ylim = c(0, 30))
# # lines(shortreperfusionOXPHOS$t, dpH, cex = 0.1, col = "red")
# plot(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$dPsi+dpH, col = "red",
#      xlab = "Time (s)", ylab = "Proton Motive Force (mV)", ylim = c(150, 200),
#      cex = 0.1)
# lines(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$dPsi+dpH, col = "red")
# 
# 
# plot(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$ATP_c, cex = 0.1,
#      col = "red")


# dpH <- -0.059*log(shortreperfusionOXPHOS2$H_x/shortreperfusionOXPHOS2$H_i)*1000
# 
# plot(shortreperfusionOXPHOS2$t, shortreperfusionOXPHOS2$NADH_x/9.4e-3,
#      xlab = "Time (s)", ylab = "NADH Proportion of Nicotinamide adenine Pool",
#      cex = 0.1, ylim = c(0, 1), col = "red")
# lines(shortreperfusionOXPHOS2$t, shortreperfusionOXPHOS2$NADH_x/9.4e-3, col = "red")
# plot(shortreperfusionOXPHOS2$t, shortreperfusionOXPHOS2$QH2_x/6.49e-3,
#      xlab = "Time (s)", ylab = "Proportion of Coenzyme Q in a Reduced State",
#      cex = 0.1, ylim = c(0, 1), col = "red")
# lines(shortreperfusionOXPHOS2$t, shortreperfusionOXPHOS2$QH2_x/6.49e-3, col = "red")
# plot(shortreperfusionOXPHOS2$t, shortreperfusionOXPHOS2$dPsi,
#      xlab = "Time (s)", ylab = "Electrical Potential Gradient (mV)",
#      cex = 0.1, col = "red")
# lines(shortreperfusionOXPHOS2$t, shortreperfusionOXPHOS2$dPsi, col = "red")
# plot(shortreperfusionOXPHOS2$t, dpH, cex = 0.1, xlab = "Time (s)", 
#      ylab= "pH Difference (mV)", col = "red", ylim = c(0, 30))
# lines(shortreperfusionOXPHOS2$t, dpH, cex = 0.1, col = "red")

# shortReperfusionPool1 <- reperfusionPool1[reperfusionPool1$t < 250, ]
# par(mfrow = c(1, 2))
# plot(shortReperfusionPool1$t, shortReperfusionPool1$NADH_x/4.7e-3, col = "red",
#      xlab = "Time (s)", ylab = "NADH fraction of NADH/NAD+ Pool", cex = 0.1,
#      ylim = c(0, 1))
# lines(shortReperfusionPool1$t, shortReperfusionPool1$NADH_x/4.7e-3, col = "red",
#      cex = 0.1)
# plot(shortReperfusionPool1$t, shortReperfusionPool1$QH2_x/6.49e-3, col = "red",
#      xlab = "Time (s)", ylab = "Reduced Proportion of Coenzyme Q Pool", cex = 0.1,
#      ylim = c(0, 1))
# lines(shortReperfusionPool1$t, shortReperfusionPool1$QH2_x/6.49e-3, col = "red",
#       cex = 0.1)
# par(mfrow = c(1, 1))
# # plot(shortReperfusionPool1$t, shortReperfusionPool1$Cred_i/4.39e-3, col = "red",
# #      xlab = "Time (s)", ylab = "Reduced Proportion of cytochrome C Pool", cex = 0.1,
# #      ylim = c(0, 1))
# # lines(shortReperfusionPool1$t, shortReperfusionPool1$Cred_i/4.39e-3, col = "red",
# #       cex = 0.1)

# shortReperfusionPoolOXPHOS1 <-
#   reperfusionPoolOXPHOS1[reperfusionPoolOXPHOS1$t < 250, ]
# dpH <- -0.059*log(shortReperfusionPoolOXPHOS1$H_x/shortReperfusionPoolOXPHOS1$H_i)*1000
# 
# par(mfrow = c(1, 2))
# plot(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$NADH_x/4.85e-3,
#      col = "red",
#      xlab = "Time (s)", ylab = "NADH fraction of NADH/NAD+ Pool", cex = 0.1,
#      ylim = c(0, 1))
# lines(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$NADH_x/4.85e-3,
#       col = "red",
#       cex = 0.1)
# plot(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$QH2_x/6.49e-3,
#      col = "red",
#      xlab = "Time (s)", ylab = "Reduced Proportion of Coenzyme Q Pool", cex = 0.1,
#      ylim = c(0, 1))
# lines(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$QH2_x/6.49e-3,
#       col = "red",
#       cex = 0.1)
# plot(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$dPsi,
#      col = "red", xlab = "Time (s)", ylab = "Electrical Potential Gradient (mV)",
#      cex = 0.1)
# lines(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$dPsi, col = "red")

# par(mfrow = c(1, 1))
# plot(shortReperfusionPoolOXPHOS1$t, 
#      shortReperfusionPoolOXPHOS1$dPsi+dpH, col = "red", xlab = "Time (s)", 
#      ylab = "Proton Motive Force (mV)", cex = 0.1)
# lines(shortReperfusionPoolOXPHOS1$t, 
#       shortReperfusionPoolOXPHOS1$dPsi+dpH, col = "red")
# # plot(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$Cred_i/4.39e-3, col = "red",
# #      xlab = "Time (s)", ylab = "Reduced Proportion of cytochrome C Pool", cex = 0.1,
# #      ylim = c(0, 1))
# # lines(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$Cred_i/4.39e-3, col = "red",
# #       cex = 0.1)


# shortReperfusionOXPHOS12 <- reperfusionOXPHOS12[
#   reperfusionOXPHOS12$t < 250, ]
# 
# dpH <- -0.059*log(shortReperfusionOXPHOS12$H_x/shortReperfusionOXPHOS12$H_i)*1000
# 
# plot(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$NADH_x/9.7e-3,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#      ylab = "Proportion of NADH/NAD+ Pool in Reduced State")
# lines(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$NADH_x/9.7e-3,
#       col = "red")
# plot(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$QH2_x/6.49e-3,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#      ylab = "Proportion of Coenzyme Q Pool in Reduced State")
# lines(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$QH2_x/6.49e-3,
#       col = "red")
# plot(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$Cred_i/4.39e-3,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#      ylab = "Proportion of Cytochrome C Pool in Reduced State")
# lines(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$Cred_i/4.39e-3,
#       col = "red")
# plot(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$dPsi+dpH,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(150, 200),
#      ylab = "Proton Motive Force (mV)")
# lines(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$dPsi+dpH,
#       col = "red")


shortReperfusionPoolOXPHOS12 <- reperfusionPoolOXPHOS12[
  reperfusionPoolOXPHOS12$t < 200, ]
dpH <- -0.059*log(shortReperfusionPoolOXPHOS12$H_x/shortReperfusionPoolOXPHOS12$H_i)*1000

## Figure 4.5

pdf("../dataVis/ischemiaRealistic.pdf", width = 15)
par(mfrow = c(1, 4), cex.axis = 2, cex.lab = 2.5, mar = c(5.1, 5, 4.1, 3))
plot(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$NADH_x/4.85e-3,
     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
     ylab = "Proportion of NADH/NAD+ Pool in Reduced State")
lines(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$NADH_x/4.85e-3,
      col = "red")
plot(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$QH2_x/6.49e-3,
     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
     ylab = "Proportion of Coenzyme Q Pool in Reduced State")
lines(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$QH2_x/6.49e-3,
      col = "red")
# plot(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$Cred_i/4.39e-3,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#      ylab = "Proportion of Cytochrome C Pool in Reduced State")
# lines(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$Cred_i/4.39e-3,
#       col = "red")
plot(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$dPsi+dpH,
     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(100, 250),
     ylab = "Proton Motive Force (mV)")
lines(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$dPsi+dpH,
      col = "red")
plot(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$SUC_x*1000, cex = 0.1,
     col = "red", xlab = "Time (s)", ylab = "Matrix Succinate (mM)")
lines(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$SUC_x*1000,
      col = "red")
par(mfrow = c(1, 1), cex.axis = 1, cex.lab = 1)
dev.off()

# shortReperfusionOXPHOS12 <- reperfusionOXPHOS12[
#   reperfusionOXPHOS12$t < 250, ]
# 
# dpH <- -0.059*log(shortReperfusionOXPHOS12$H_x/shortReperfusionOXPHOS12$H_i)*1000
# 
# plot(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$NADH_x/9.7e-3,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#      ylab = "Proportion of NADH/NAD+ Pool in Reduced State")
# lines(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$NADH_x/9.7e-3,
#       col = "red")
# plot(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$QH2_x/6.49e-3,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#      ylab = "Proportion of Coenzyme Q Pool in Reduced State")
# lines(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$QH2_x/6.49e-3,
#       col = "red")
# plot(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$Cred_i/4.39e-3,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#      ylab = "Proportion of Cytochrome C Pool in Reduced State")
# lines(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$Cred_i/4.39e-3,
#       col = "red")
# plot(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$dPsi+dpH,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(150, 200),
#      ylab = "Proton Motive Force (mV)")
# lines(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$dPsi+dpH,
#       col = "red")


shortReperfusionPool12 <- reperfusionPool12[
  reperfusionPool12$t < 200, ]
dpH <- -0.059*log(shortReperfusionPool12$H_x/shortReperfusionPool12$H_i)*1000


pdf("../dataVis/ischemiaReperfusionPoolRedox.pdf", width = 15)
par(mfrow = c(1, 4), cex.axis = 2, cex.lab = 2.5, mar = c(5.1, 5, 4.1, 3))
plot(shortReperfusionPool12$t, shortReperfusionPool12$NADH_x/4.85e-3,
     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
     ylab = "Proportion of NADH/NAD+ Pool in Reduced State")
lines(shortReperfusionPool12$t, shortReperfusionPool12$NADH_x/4.85e-3,
      col = "red")
plot(shortReperfusionPool12$t, shortReperfusionPool12$QH2_x/6.49e-3,
     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
     ylab = "Proportion of Coenzyme Q Pool in Reduced State")
lines(shortReperfusionPool12$t, shortReperfusionPool12$QH2_x/6.49e-3,
      col = "red")
# plot(shortReperfusionPool12$t, shortReperfusionPool12$Cred_i/4.39e-3,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#      ylab = "Proportion of Cytochrome C Pool in Reduced State")
# lines(shortReperfusionPool12$t, shortReperfusionPool12$Cred_i/4.39e-3,
#       col = "red")
#pdf("../dataVis/pmfIschemiaReperfusionPool.pdf")
#par(cex.axis = 2, cex.lab = 2.5, mar = c(5.1, 5, 4.1, 2.1))
plot(shortReperfusionPool12$t, shortReperfusionPool12$dPsi+dpH,
     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(100, 260),
     ylab = "Proton Motive Force (mV)")
lines(shortReperfusionPool12$t, shortReperfusionPool12$dPsi+dpH,
      col = "red")
plot(shortReperfusionPool12$t, shortReperfusionPool12$SUC_x*1000, cex = 0.1,
     col = "red", xlab = "Time (s)", ylab = "Matrix Succinate (mM)")
lines(shortReperfusionPool12$t, shortReperfusionPool12$SUC_x*1000,
      col = "red")
par(mfrow = c(1, 1))
dev.off()
