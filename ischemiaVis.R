if (!grepl("modelScripts", getwd())) {
   setwd("modelScripts")
}

# reperfusionOXPHOS <- feather::read_feather("../results/resultsReperfusionOXPHOS1.feather")
# reperfusionOXPHOS2 <- feather::read_feather("../results/resultsReperfusionOXPHOS2.feather")
# reperfusionOXPHOS12 <- feather::read_feather("../results/resultsReperfusionOXPHOS12.feather")

# reperfusionPool1 <- feather::read_feather(
#    "../results/resultsReperfusionPool1.feather")
reperfusionPool12 <- data.table::fread(
   "../results/resultsReperfusionPool12.csv")

# reperfusionPoolOXPHOS1 <- feather::read_feather(
#    "../results/resultsReperfusionOXPHOSPool1.feather")
reperfusionPoolOXPHOS12 <- data.table::fread(
   "../results/resultsReperfusionOXPHOSPool12.csv")

# reperfusion1 <- read.csv("../results/resultsReperfusion1.csv")
# shortreperfusion <- reperfusion1[reperfusion1$t < 120, ]

# plot(shortreperfusion$t, shortreperfusion$dPsi, cex = 0.1,
#      col = "red", ylim = c(150, 180), ylab = "dPsi (mV)",
#      xlab = "Time (s)")
# lines(shortreperfusion$t, shortreperfusion$dPsi, cex = 0.1,
#       col = "red")

# reperfusion2 <- read.csv("../results/resultsReperfusion2.csv")
# shortreperfusion <- reperfusion2[reperfusion2$t < 120, ]

# plot(shortreperfusion$t, shortreperfusion$dPsi, cex = 0.1,
#      col = "red", ylim = c(150, 180), ylab = "dPsi (mV)",
#      xlab = "Time (s)")
# lines(shortreperfusion$t, shortreperfusion$dPsi, cex = 0.1,
#       col = "red")

# reperfusion12 <- read.csv("../results/resultsReperfusion12.csv")
# shortreperfusion12 <- reperfusion12[reperfusion12$t < 120, ]

# plot(shortreperfusion$t, shortreperfusion$dPsi, cex = 0.1,
#      col = "red", ylim = c(150, 180), ylab = "dPsi (mV)",
#      xlab = "Time (s)")
# lines(shortreperfusion$t, shortreperfusion$dPsi, cex = 0.1,
#       col = "red")

# shortreperfusion1 <- reperfusion1[reperfusion1$t < 60, ]
# shortreperfusion2 <- reperfusion2[reperfusion2$t < 60, ]
# shortreperfusion12 <- reperfusion12[reperfusion12$t < 60, ]
# 
# pdf("../dataVis/cytCIR.pdf")
# par(mfrow = c(1, 3), cex.lab = 1.5)
# 
# plot(shortreperfusion1$t, shortreperfusion1$Cred_i/1.956e-3, cex = 0.1,
#      col = "red", ylim = c(0, 1), ylab = "Cytochrome C Reduced Proportion",
#      xlab = "Time (s)")
# lines(shortreperfusion1$t, shortreperfusion1$Cred_i/1.956e-3, cex = 0.1,
#       col = "red")
# 
# plot(shortreperfusion2$t, shortreperfusion2$Cred_i/1.956e-3, cex = 0.1,
#      col = "red", ylim = c(0, 1), ylab = "",
#      xlab = "Time (s)")
# lines(shortreperfusion2$t, shortreperfusion2$Cred_i/1.956e-3, cex = 0.1,
#       col = "red")
# 
# plot(shortreperfusion12$t, shortreperfusion12$Cred_i/1.956e-3, cex = 0.1,
#      col = "red", ylim = c(0, 1), ylab = "",
#      xlab = "Time (s)")
# lines(shortreperfusion12$t, shortreperfusion12$Cred_i/1.956e-3, cex = 0.1,
#       col = "red")
# 
# dev.off()
# 
# pdf("../dataVis/coqIR.pdf")
# par(mfrow = c(1, 3), cex.lab = 1.5, cex.axis = 1.5)
# plot(shortreperfusion1$t, shortreperfusion1$QH2_x/2.148e-3, cex = 0.1,
#      col = "red", ylim = c(0, 1), ylab = "Coenzyme Q Reduced Proportion",
#      xlab = "Time (s)")
# lines(shortreperfusion1$t, shortreperfusion1$QH2_x/2.148e-3, cex = 0.1,
#       col = "red")
# 
# plot(shortreperfusion2$t, shortreperfusion2$QH2_x/2.148e-3, cex = 0.1,
#      col = "red", ylim = c(0, 1), ylab = "",
#      xlab = "Time (s)")
# lines(shortreperfusion2$t, shortreperfusion2$QH2_x/2.148e-3, cex = 0.1,
#       col = "red")
# 
# plot(shortreperfusion12$t, shortreperfusion12$QH2_x/2.148e-3, cex = 0.1,
#      col = "red", ylim = c(0, 1), ylab = "",
#      xlab = "Time (s)")
# lines(shortreperfusion12$t, shortreperfusion12$QH2_x/2.148e-3, cex = 0.1,
#       col = "red")
# dev.off()
# 
# pdf("../dataVis/nadhIR.pdf")
# par(mfrow = c(1, 3), cex.lab = 1.5, cex.axis = 1.5)
# plot(shortreperfusion1$t, shortreperfusion1$NADH_x/0.824e-3, cex = 0.1,
#      col = "red", ylim = c(0, 1), ylab = "NADH Reduced Proportion",
#      xlab = "Time (s)")
# lines(shortreperfusion1$t, shortreperfusion1$NADH_x/0.824e-3, cex = 0.1,
#       col = "red")
# 
# plot(shortreperfusion2$t, shortreperfusion2$NADH_x/0.824e-3, cex = 0.1,
#      col = "red", ylim = c(0, 1), ylab = "",
#      xlab = "Time (s)")
# lines(shortreperfusion2$t, shortreperfusion2$NADH_x/0.824e-3, cex = 0.1,
#       col = "red")
# 
# plot(shortreperfusion12$t, shortreperfusion12$NADH_x/0.824e-3, cex = 0.1,
#      col = "red", ylim = c(0, 1), ylab = "",
#      xlab = "Time (s)")
# lines(shortreperfusion12$t, shortreperfusion12$NADH_x/0.824e-3, cex = 0.1,
#       col = "red")
# dev.off()
# 
# plot(shortreperfusion1$t, shortreperfusion1$H_x, cex = 0.1, col = "red",
#      ylab = "", xlab = "Time (s)", ylim = c(3.8e-8, 5.5e-8))
# 
# plot(shortreperfusion2$t, shortreperfusion2$H_x, cex = 0.1, col = "red",
#      ylab = "", xlab = "Time (s)", ylim = c(3.8e-8, 5.5e-8))
# 
# plot(shortreperfusion12$t, shortreperfusion12$H_x, cex = 0.1, col = "red",
#      ylab = "", xlab = "Time (s)", ylim = c(3.8e-8, 5.5e-8))
# 
# par(mfrow = c(1, 1))
# 
# shortreperfusionOXPHOS <- reperfusionOXPHOS[reperfusionOXPHOS$t < 20, ]
# shortreperfusionOXPHOS2 <- reperfusionOXPHOS2[reperfusionOXPHOS2$t < 20, ]
# 
# dpH <- -0.059*log(shortreperfusionOXPHOS$H_x/shortreperfusionOXPHOS$H_i)*1000
# 
# par(mfrow = c(1, 2))
# plot(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$NADH_x/8.24e-4,
#      xlab = "Time (s)", ylab = "NADH Proportion of Nicotinamide adenine Pool",
#      cex = 0.1, ylim = c(0, 1), col = "red")
# lines(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$NADH_x/8.24e-4, col = "red")
# plot(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$QH2_x/2.148e-3,
#      xlab = "Time (s)", ylab = "Proportion of Coenzyme Q in a Reduced State",
#      cex = 0.1, ylim = c(0, 1), col = "red")
# lines(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$QH2_x/2.148e-3, col = "red")
# par(mfrow = c(1, 1))
# #plot(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$dPsi + dpH,
# #     xlab = "Time (s)", ylab = "Electrical Potential Gradient (mV)",
# #     cex = 0.1, col = "red", ylim = c(140, 180))
# #lines(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$dPsi + dpH, col = "red")
# # plot(shortreperfusionOXPHOS$t, dpH, cex = 0.1, xlab = "Time (s)", 
# #      ylab= "pH Difference (mV)", col = "red", ylim = c(0, 30))
# # lines(shortreperfusionOXPHOS$t, dpH, cex = 0.1, col = "red")
# plot(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$dPsi+dpH, col = "red",
#      xlab = "Time (s)", ylab = "Proton Motive Force (mV)", ylim = c(150, 170),
#      cex = 0.1)
# lines(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$dPsi+dpH, col = "red")
# 
# 
# plot(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$ATP_c, cex = 0.1,
#      col = "red")
# 
# 
# dpH <- -0.059*log(shortreperfusionOXPHOS2$H_x/shortreperfusionOXPHOS2$H_i)*1000
#  
# plot(shortreperfusionOXPHOS2$t, shortreperfusionOXPHOS2$NADH_x/0.824e-3,
#      xlab = "Time (s)", ylab = "NADH Proportion of Nicotinamide adenine Pool",
#      cex = 0.1, ylim = c(0, 1), col = "red")
# lines(shortreperfusionOXPHOS2$t, shortreperfusionOXPHOS2$NADH_x/0.824e-3, col = "red")
# # plot(shortreperfusionOXPHOS2$t, shortreperfusionOXPHOS2$QH2_x/2.148e-3,
# #      xlab = "Time (s)", ylab = "Proportion of Coenzyme Q in a Reduced State",
# #      cex = 0.1, ylim = c(0, 1), col = "red")
# # lines(shortreperfusionOXPHOS2$t, shortreperfusionOXPHOS2$QH2_x/2.148e-3, col = "red")
# # plot(shortreperfusionOXPHOS2$t, shortreperfusionOXPHOS2$dPsi,
# #      xlab = "Time (s)", ylab = "Electrical Potential Gradient (mV)",
# #      cex = 0.1, col = "red")
# # lines(shortreperfusionOXPHOS2$t, shortreperfusionOXPHOS2$dPsi, col = "red")
# # plot(shortreperfusionOXPHOS2$t, dpH, cex = 0.1, xlab = "Time (s)", 
# #      ylab= "pH Difference (mV)", col = "red", ylim = c(0, 30))
# # lines(shortreperfusionOXPHOS2$t, dpH, cex = 0.1, col = "red")
# 
# shortReperfusionPool1 <- reperfusionPool1[reperfusionPool1$t < 100, ]
# par(mfrow = c(1, 2))
# plot(shortReperfusionPool1$t, shortReperfusionPool1$NADH_x/4.7e-3, col = "red",
#      xlab = "Time (s)", ylab = "NADH fraction of NADH/NAD+ Pool", cex = 0.1,
#      ylim = c(0, 1))
# lines(shortReperfusionPool1$t, shortReperfusionPool1$NADH_x/4.7e-3, col = "red",
#       cex = 0.1)
# plot(shortReperfusionPool1$t, shortReperfusionPool1$QH2_x/2.148e-3, col = "red",
#      xlab = "Time (s)", ylab = "Reduced Proportion of Coenzyme Q Pool", cex = 0.1,
#      ylim = c(0, 1))
# lines(shortReperfusionPool1$t, shortReperfusionPool1$QH2_x/2.148e-3, col = "red",
#       cex = 0.1)
# par(mfrow = c(1, 1))
# # plot(shortReperfusionPool1$t, shortReperfusionPool1$Cred_i/1.956e-3, col = "red",
# #      xlab = "Time (s)", ylab = "Reduced Proportion of cytochrome C Pool", cex = 0.1,
# #      ylim = c(0, 1))
# # lines(shortReperfusionPool1$t, shortReperfusionPool1$Cred_i/1.956e-3, col = "red",
# #       cex = 0.1)
# 
# shortReperfusionPoolOXPHOS1 <- 
#    reperfusionPoolOXPHOS1[reperfusionPoolOXPHOS1$t < 200, ]
# dpH <- -0.059*log(shortReperfusionPoolOXPHOS1$H_x/shortReperfusionPoolOXPHOS1$H_i)*1000
# 
# par(mfrow = c(1, 2))
# plot(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$NADH_x/8.24e-4, 
#      col = "red",
#      xlab = "Time (s)", ylab = "NADH fraction of NADH/NAD+ Pool", cex = 0.1,
#      ylim = c(0, 1))
# lines(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$NADH_x/8.2e-4, 
#       col = "red",
#       cex = 0.1)
# plot(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$QH2_x/2.148e-3, 
#      col = "red",
#      xlab = "Time (s)", ylab = "Reduced Proportion of Coenzyme Q Pool", cex = 0.1,
#      ylim = c(0, 1))
# lines(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$QH2_x/2.148e-3, 
#       col = "red",
#       cex = 0.1)
# # plot(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$dPsi,
# #      col = "red", xlab = "Time (s)", ylab = "Electrical Potential Gradient (mV)",
# #      cex = 0.1)
# # lines(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$dPsi, col = "red")
# 
# par(mfrow = c(1, 1))
# plot(shortReperfusionPoolOXPHOS1$t, 
#      shortReperfusionPoolOXPHOS1$dPsi+dpH, col = "red", xlab = "Time (s)", 
#      ylab = "Proton Motive Force (mV)", cex = 0.1)
# lines(shortReperfusionPoolOXPHOS1$t, 
#       shortReperfusionPoolOXPHOS1$dPsi+dpH, col = "red")
# # plot(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$Cred_i/1.956e-3, col = "red",
# #      xlab = "Time (s)", ylab = "Reduced Proportion of cytochrome C Pool", cex = 0.1,
# #      ylim = c(0, 1))
# # lines(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$Cred_i/1.956e-3, col = "red",
# #       cex = 0.1)
# 
# 
# shortReperfusionOXPHOS12 <- reperfusionOXPHOS12[
#    reperfusionOXPHOS12$t < 20, ]
# 
# dpH <- -0.059*log(shortReperfusionOXPHOS12$H_x/shortReperfusionOXPHOS12$H_i)*1000
# 
# plot(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$NADH_x/0.824e-3,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#      ylab = "Proportion of NADH/NAD+ Pool in Reduced State")
# lines(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$NADH_x/0.824e-3,
#       col = "red")
# plot(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$QH2_x/2.148e-3,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#      ylab = "Proportion of Coenzyme Q Pool in Reduced State")
# lines(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$QH2_x/2.148e-3,
#       col = "red")
# plot(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$Cred_i/1.956e-3,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#      ylab = "Proportion of Cytochrome C Pool in Reduced State")
# lines(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$Cred_i/1.956e-3,
#       col = "red")
# plot(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$dPsi+dpH,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(150, 200),
#      ylab = "Proton Motive Force (mV)")
# lines(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$dPsi+dpH,
#       col = "red")


shortReperfusionPoolOXPHOS12 <- reperfusionPoolOXPHOS12[
   reperfusionPoolOXPHOS12$t < 20, ]
dpH <- -0.059*log(shortReperfusionPoolOXPHOS12$H_x/shortReperfusionPoolOXPHOS12$H_i)*1000

plot(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$NADH_x/0.824e-3,
     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
     ylab = "Proportion of NADH/NAD+ Pool in Reduced State")
lines(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$NADH_x/0.824e-3,
      col = "red")

plot(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$QH2_x/2.148e-3,
     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
     ylab = "Proportion of Coenzyme Q Pool in Reduced State")
lines(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$QH2_x/2.148e-3,
      col = "red")
plot(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$Cred_i/1.956e-3,
     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
     ylab = "Proportion of Cytochrome C Pool in Reduced State")
lines(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$Cred_i/1.956e-3,
      col = "red")
plot(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$dPsi+dpH,
     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(150, 210),
     ylab = "Proton Motive Force (mV)")
lines(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$dPsi+dpH,
      col = "red")

# 
# shortReperfusionOXPHOS12 <- reperfusionOXPHOS12[
#    reperfusionOXPHOS12$t < 150, ]
# 
# dpH <- -0.059*log(shortReperfusionOXPHOS12$H_x/shortReperfusionOXPHOS12$H_i)*1000
# 
# plot(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$NADH_x/0.824e-3,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#      ylab = "Proportion of NADH/NAD+ Pool in Reduced State")
# lines(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$NADH_x/0.824e-3,
#       col = "red")
# plot(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$QH2_x/2.148e-3,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#      ylab = "Proportion of Coenzyme Q Pool in Reduced State")
# lines(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$QH2_x/2.148e-3,
#       col = "red")
# plot(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$Cred_i/1.956e-3,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#      ylab = "Proportion of Cytochrome C Pool in Reduced State")
# lines(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$Cred_i/1.956e-3,
#       col = "red")
# plot(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$dPsi+dpH,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(150, 200),
#      ylab = "Proton Motive Force (mV)")
# lines(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$dPsi+dpH,
#       col = "red")


shortReperfusionPool12 <- reperfusionPool12[
   reperfusionPool12$t < 150, ]
dpH <- -0.059*log(shortReperfusionPool12$H_x/shortReperfusionPool12$H_i)*1000

pdf("../dataVis/ischemiaReperfusionPoolPT.pdf")
par(mfrow = c(1, 2))
## Important case
plot(shortReperfusionPool12$t, shortReperfusionPool12$NADH_x/0.824e-3,
     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
     ylab = "Proportion of NADH/NAD+ Pool in Reduced State")
lines(shortReperfusionPool12$t, shortReperfusionPool12$NADH_x/0.824e-3,
      col = "red")
# plot(shortReperfusionPool12$t, shortReperfusionPool12$QH2_x/2.148e-3,
#     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#     ylab = "Proportion of Coenzyme Q Pool in Reduced State")
# lines(shortReperfusionPool12$t, shortReperfusionPool12$QH2_x/2.148e-3,
#      col = "red")
# plot(shortReperfusionPool12$t, shortReperfusionPool12$Cred_i/1.956e-3,
#     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#     ylab = "Proportion of Cytochrome C Pool in Reduced State")
# lines(shortReperfusionPool12$t, shortReperfusionPool12$Cred_i/1.956e-3,
#      col = "red")
plot(shortReperfusionPool12$t, shortReperfusionPool12$dPsi+dpH,
     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 200),
     ylab = "Proton Motive Force (mV)")
lines(shortReperfusionPool12$t, shortReperfusionPool12$dPsi+dpH,
      col = "red")
par(mfrow = c(1, 1))
dev.off()
