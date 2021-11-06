#reperfusion1 <- read.csv("../results/resultsmTALReperfusion1.csv")
#shortreperfusion <- reperfusion1[reperfusion1$t < 120, ]

# plot(shortreperfusion$t, shortreperfusion$dPsi, cex = 0.1,
#      col = "red", ylim = c(150, 180), ylab = "dPsi (mV)",
#      xlab = "Time (s)")
# lines(shortreperfusion$t, shortreperfusion$dPsi, cex = 0.1,
#       col = "red")

#reperfusion2 <- read.csv("../results/resultsmTALReperfusion2.csv")
#shortreperfusion <- reperfusion2[reperfusion2$t < 120, ]

# plot(shortreperfusion$t, shortreperfusion$dPsi, cex = 0.1,
#      col = "red", ylim = c(150, 180), ylab = "dPsi (mV)",
#      xlab = "Time (s)")
# lines(shortreperfusion$t, shortreperfusion$dPsi, cex = 0.1,
#       col = "red")

#reperfusion12 <- read.csv("../results/resultsmTALReperfusion12.csv")
#shortreperfusion12 <- reperfusion12[reperfusion12$t < 120, ]

# plot(shortreperfusion$t, shortreperfusion$dPsi, cex = 0.1,
#      col = "red", ylim = c(150, 180), ylab = "dPsi (mV)",
#      xlab = "Time (s)")
# lines(shortreperfusion$t, shortreperfusion$dPsi, cex = 0.1,
#       col = "red")

#shortreperfusion1 <- reperfusion1[reperfusion1$t < 10, ]
#shortreperfusion2 <- reperfusion2[reperfusion2$t < 10, ]
#shortreperfusion12 <- reperfusion12[reperfusion12$t < 10, ]

#par(mfrow = c(1, 3))

#plot(shortreperfusion1$t, shortreperfusion1$Cred_i, cex = 0.1,
#     col = "red", ylim = c(0, 0.0016), ylab = "Reduced Cytochrome C (M)",
#     xlab = "Time (s)")
#lines(shortreperfusion1$t, shortreperfusion1$Cred_i, cex = 0.1,
#      col = "red")

#plot(shortreperfusion2$t, shortreperfusion2$Cred_i, cex = 0.1,
#     col = "red", ylim = c(0, 0.0016), ylab = "",
#     xlab = "Time (s)")
#lines(shortreperfusion2$t, shortreperfusion2$Cred_i, cex = 0.1,
#      col = "red")

#plot(shortreperfusion12$t, shortreperfusion12$Cred_i, cex = 0.1,
#     col = "red", ylim = c(0, 0.0016), ylab = "",
#     xlab = "Time (s)")
#lines(shortreperfusion12$t, shortreperfusion12$Cred_i, cex = 0.1,
#      col = "red")

#plot(shortreperfusion1$t, shortreperfusion1$QH2_x, cex = 0.1,
#     col = "red", ylim = c(0, 0.0014), ylab = "Reduced Coenzyme Q (M)",
#     xlab = "Time (s)")
#lines(shortreperfusion1$t, shortreperfusion1$QH2_x, cex = 0.1,
#      col = "red")

#plot(shortreperfusion2$t, shortreperfusion2$QH2_x, cex = 0.1,
#     col = "red", ylim = c(0, 0.0014), ylab = "",
#     xlab = "Time (s)")
#lines(shortreperfusion2$t, shortreperfusion2$QH2_x, cex = 0.1,
#      col = "red")

#plot(shortreperfusion12$t, shortreperfusion12$QH2_x, cex = 0.1,
#     col = "red", ylim = c(0, 0.0014), ylab = "",
#     xlab = "Time (s)")
#lines(shortreperfusion12$t, shortreperfusion12$QH2_x, cex = 0.1,
#      col = "red")

#plot(shortreperfusion1$t, shortreperfusion1$NADH_x, cex = 0.1,
#     col = "red", ylim = c(0, 0.0001), ylab = "NADH (M)",
#     xlab = "Time (s)")
#lines(shortreperfusion1$t, shortreperfusion1$NADH_x, cex = 0.1,
#      col = "red")

#plot(shortreperfusion2$t, shortreperfusion2$NADH_x, cex = 0.1,
#     col = "red", ylab = "", ylim = c(0, 0.0001),
#     xlab = "Time (s)")
#lines(shortreperfusion2$t, shortreperfusion2$NADH_x, cex = 0.1,
#      col = "red")

#plot(shortreperfusion12$t, shortreperfusion12$NADH_x, cex = 0.1,
#     col = "red", ylab = "", ylim = c(0, 0.0001),
#     xlab = "Time (s)")
#lines(shortreperfusion12$t, shortreperfusion12$NADH_x, cex = 0.1,
#      col = "red")

par(mfrow = c(1, 1))

#reperfusionOXPHOS <- feather::read_feather("../results/resultsReperfusionOXPHOS1mTAL.feather")
#reperfusionOXPHOS2 <- feather::read_feather("../results/resultsReperfusionOXPHOS2mTAL.feather")
#reperfusionOXPHOS12 <- feather::read_feather("../results/resultsReperfusionOXPHOS12mTAL.feather")

#reperfusionPool1 <- feather::read_feather(
#   "../results/resultsReperfusionPool1mTAL.feather")
reperfusionPool12 <- data.table::fread(
   "../results/resultsReperfusionPool12mTAL.csv")

#reperfusionPoolOXPHOS1 <- feather::read_feather(
#   "../results/resultsReperfusionOXPHOSPool1mTAL.feather")
reperfusionPoolOXPHOS12 <- data.table::fread(
   "../results/resultsReperfusionOXPHOSPool12mTAL.csv")


#shortreperfusionOXPHOS <- reperfusionOXPHOS[reperfusionOXPHOS$t < 20, ]
#shortreperfusionOXPHOS2 <- reperfusionOXPHOS2[reperfusionOXPHOS2$t < 20, ]

#dpH <- -0.059*log(shortreperfusionOXPHOS$H_x/shortreperfusionOXPHOS$H_i)*1000

#par(mfrow = c(1, 2))
#plot(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$NADH_x/8.24e-5,
#     xlab = "Time (s)", ylab = "NADH Proportion of Nicotinamide adenine Pool",
#     cex = 0.1, ylim = c(0, 1), col = "red")
#lines(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$NADH_x/8.24e-5, col = "red")
#plot(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$QH2_x/2.148e-3,
#     xlab = "Time (s)", ylab = "Proportion of Coenzyme Q in a Reduced State",
#     cex = 0.1, ylim = c(0, 1), col = "red")
#lines(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$QH2_x/2.148e-3, col = "red")
#par(mfrow = c(1, 1))
# plot(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$dPsi,
#      xlab = "Time (s)", ylab = "Electrical Potential Gradient (mV)",
#      cex = 0.1, col = "red", ylim = c(140, 180))
# lines(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$dPsi, col = "red")
# plot(shortreperfusionOXPHOS$t, dpH, cex = 0.1, xlab = "Time (s)", 
#      ylab= "pH Difference (mV)", col = "red", ylim = c(0, 30))
# lines(shortreperfusionOXPHOS$t, dpH, cex = 0.1, col = "red")
#plot(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$dPsi+dpH, col = "red",
#     xlab = "Time (s)", ylab = "Proton Motive Force (mV)", ylim = c(150, 170),
#     cex = 0.1)
#lines(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$dPsi+dpH, col = "red")


#plot(shortreperfusionOXPHOS$t, shortreperfusionOXPHOS$ATP_c, cex = 0.1,
#     col = "red")


# dpH <- -0.059*log(shortreperfusionOXPHOS2$H_x/shortreperfusionOXPHOS2$H_i)*1000
# 
# plot(shortreperfusionOXPHOS2$t, shortreperfusionOXPHOS2$NADH_x/9.4e-3,
#      xlab = "Time (s)", ylab = "NADH Proportion of Nicotinamide adenine Pool",
#      cex = 0.1, ylim = c(0, 1), col = "red")
# lines(shortreperfusionOXPHOS2$t, shortreperfusionOXPHOS2$NADH_x/9.4e-3, col = "red")
# plot(shortreperfusionOXPHOS2$t, shortreperfusionOXPHOS2$QH2_x/2.148e-3,
#      xlab = "Time (s)", ylab = "Proportion of Coenzyme Q in a Reduced State",
#      cex = 0.1, ylim = c(0, 1), col = "red")
# lines(shortreperfusionOXPHOS2$t, shortreperfusionOXPHOS2$QH2_x/2.148e-3, col = "red")
# plot(shortreperfusionOXPHOS2$t, shortreperfusionOXPHOS2$dPsi,
#      xlab = "Time (s)", ylab = "Electrical Potential Gradient (mV)",
#      cex = 0.1, col = "red")
# lines(shortreperfusionOXPHOS2$t, shortreperfusionOXPHOS2$dPsi, col = "red")
# plot(shortreperfusionOXPHOS2$t, dpH, cex = 0.1, xlab = "Time (s)", 
#      ylab= "pH Difference (mV)", col = "red", ylim = c(0, 30))
# lines(shortreperfusionOXPHOS2$t, dpH, cex = 0.1, col = "red")

#shortReperfusionPool1 <- reperfusionPool1[reperfusionPool1$t < 100, ]
#par(mfrow = c(1, 2))
#plot(shortReperfusionPool1$t, shortReperfusionPool1$NADH_x/0.824e-4, col = "red",
#     xlab = "Time (s)", ylab = "NADH fraction of NADH/NAD+ Pool", cex = 0.1,
#     ylim = c(0, 1))
#lines(shortReperfusionPool1$t, shortReperfusionPool1$NADH_x/0.824e-4, col = "red",
#      cex = 0.1)
#plot(shortReperfusionPool1$t, shortReperfusionPool1$QH2_x/2.148e-3, col = "red",
#     xlab = "Time (s)", ylab = "Reduced Proportion of Coenzyme Q Pool", cex = 0.1,
#     ylim = c(0, 1))
#lines(shortReperfusionPool1$t, shortReperfusionPool1$QH2_x/2.148e-3, col = "red",
#      cex = 0.1)
#par(mfrow = c(1, 1))
# plot(shortReperfusionPool1$t, shortReperfusionPool1$Cred_i/1.956e-3, col = "red",
#      xlab = "Time (s)", ylab = "Reduced Proportion of cytochrome C Pool", cex = 0.1,
#      ylim = c(0, 1))
# lines(shortReperfusionPool1$t, shortReperfusionPool1$Cred_i/1.956e-3, col = "red",
#       cex = 0.1)

#shortReperfusionPoolOXPHOS1 <- 
#   reperfusionPoolOXPHOS1[reperfusionPoolOXPHOS1$t < 200, ]
#dpH <- -0.059*log(shortReperfusionPoolOXPHOS1$H_x/shortReperfusionPoolOXPHOS1$H_i)*1000

## Important
#par(mfrow = c(1, 2))
#plot(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$NADH_x/8.24e-5, 
#     col = "red",
#     xlab = "Time (s)", ylab = "NADH fraction of NADH/NAD+ Pool", cex = 0.1,
#     ylim = c(0, 1))
#lines(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$NADH_x/8.24e-5, 
#      col = "red",
#      cex = 0.1)
#plot(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$QH2_x/2.148e-3, 
#     col = "red",
#     xlab = "Time (s)", ylab = "Reduced Proportion of Coenzyme Q Pool", cex = 0.1,
#     ylim = c(0, 1))
#lines(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$QH2_x/2.148e-3, 
#      col = "red",
#      cex = 0.1)
# plot(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$dPsi + dpH,
#      col = "red", xlab = "Time (s)", ylab = "Electrical Potential Gradient (mV)",
#      cex = 0.1)
# lines(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$dPsi + dpH, col = "red")

#par(mfrow = c(1, 1))
#plot(shortReperfusionPoolOXPHOS1$t, 
#     shortReperfusionPoolOXPHOS1$dPsi+dpH, col = "red", xlab = "Time (s)", 
#     ylab = "Proton Motive Force (mV)", cex = 0.1)
#lines(shortReperfusionPoolOXPHOS1$t, 
#      shortReperfusionPoolOXPHOS1$dPsi+dpH, col = "red")
# plot(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$Cred_i/1.956e-3, col = "red",
#      xlab = "Time (s)", ylab = "Reduced Proportion of cytochrome C Pool", cex = 0.1,
#      ylim = c(0, 1))
# lines(shortReperfusionPoolOXPHOS1$t, shortReperfusionPoolOXPHOS1$Cred_i/1.956e-3, col = "red",
#       cex = 0.1)

# 
# shortReperfusionOXPHOS12 <- reperfusionOXPHOS12[
#    reperfusionOXPHOS12$t < 20, ]
# 
# dpH <- -0.059*log(shortReperfusionOXPHOS12$H_x/shortReperfusionOXPHOS12$H_i)*1000
# 
# plot(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$NADH_x/0.824e-4,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#      ylab = "Proportion of NADH/NAD+ Pool in Reduced State")
# lines(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$NADH_x/0.824e-4,
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
   reperfusionPoolOXPHOS12$t < 50, ]
dpH <- -0.059*log(shortReperfusionPoolOXPHOS12$H_x/shortReperfusionPoolOXPHOS12$H_i)*1000

## Important case?
par(mfrow = c(1, 2))
plot(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$NADH_x/0.824e-4,
     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
     ylab = "Proportion of NADH/NAD+ Pool in Reduced State")
lines(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$NADH_x/0.824e-4,
      col = "red")
plot(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$QH2_x/2.148e-3,
     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
     ylab = "Proportion of Coenzyme Q Pool in Reduced State")
lines(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$QH2_x/2.148e-3,
      col = "red")
par(mfrow = c(1, 1))
#plot(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$Cred_i/1.956e-3,
#     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#     ylab = "Proportion of Cytochrome C Pool in Reduced State")
#lines(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$Cred_i/1.956e-3,
#      col = "red")
# plot(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$dPsi+dpH,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(150, 210),
#      ylab = "Proton Motive Force (mV)")
# lines(shortReperfusionPoolOXPHOS12$t, shortReperfusionPoolOXPHOS12$dPsi+dpH,
#       col = "red")


# shortReperfusionOXPHOS12 <- reperfusionOXPHOS12[
#    reperfusionOXPHOS12$t < 150, ]
# 
# dpH <- -0.059*log(shortReperfusionOXPHOS12$H_x/shortReperfusionOXPHOS12$H_i)*1000
# 
# plot(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$NADH_x/0.824e-4,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#      ylab = "Proportion of NADH/NAD+ Pool in Reduced State")
# lines(shortReperfusionOXPHOS12$t, shortReperfusionOXPHOS12$NADH_x/0.824e-4,
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
   reperfusionPool12$t < 200, ]
dpH <- -0.059*log(shortReperfusionPool12$H_x/shortReperfusionPool12$H_i)*1000

## Important case?
pdf("../dataVis/reperfusionPoolmTAL.pdf", width = 10)
par(mfrow = c(1, 2), cex.axis = 1.5, cex.lab = 1.5)
plot(shortReperfusionPool12$t, shortReperfusionPool12$NADH_x/0.824e-4,
     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
     ylab = "Proportion of NADH/NAD+ Pool in Reduced State")
lines(shortReperfusionPool12$t, shortReperfusionPool12$NADH_x/0.824e-4,
      col = "red")
# plot(shortReperfusionPool12$t, shortReperfusionPool12$QH2_x/2.148e-3,
#      col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#      ylab = "Proportion of Coenzyme Q Pool in Reduced State")
# lines(shortReperfusionPool12$t, shortReperfusionPool12$QH2_x/2.148e-3,
#       col = "red")
#plot(shortReperfusionPool12$t, shortReperfusionPool12$Cred_i/1.956e-3,
#     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(0, 1),
#     ylab = "Proportion of Cytochrome C Pool in Reduced State")
#lines(shortReperfusionPool12$t, shortReperfusionPool12$Cred_i/1.956e-3,
#      col = "red")
#par(mfrow = c(1, 1))

plot(shortReperfusionPool12$t, shortReperfusionPool12$dPsi+dpH,
     col = "red", cex = 0.1, xlab = "Time (s)", ylim = c(150, 200),
     ylab = "Proton Motive Force (mV)")
lines(shortReperfusionPool12$t, shortReperfusionPool12$dPsi+dpH,
      col = "red")
par(mfrow = c(1, 1))
dev.off()

