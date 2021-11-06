cyanide <- feather::read_feather("./results/resultsCyanide.feather")
baselines <- feather::read_feather("./results/resultsFe.feather")

#plot(cyanide$t, cyanide$ATP_c, cex = 0.01, xlab = "Time", ylab = "[ATP]_c")
#lines(cyanide$t, cyanide$ATP_c, cex = 0.01)

pdf("./dataVis/ATPCyanide.pdf")

plot(cyanide$t[cyanide$t < 20], cyanide$ATP_c[cyanide$t < 20], cex = 0.01,
     xlab = "Time (s)", ylab = "[ATP]_c (mM)", ylim = c(0, 0.0025), col = "red")
lines(cyanide$t[cyanide$t < 20], cyanide$ATP_c[cyanide$t < 20], cex = 0.01,
      col = "red")
points(baselines$t[baselines$t < 20], baselines$ATP_c[baselines$t < 20],
       cex = 0.01)
lines(baselines$t[baselines$t < 20], baselines$ATP_c[baselines$t < 20],
       cex = 0.01)

dev.off()

#plot(cyanide$t, cyanide$dPsi, cex = 0.01, xlab = "Time", 
#     ylab = "Mem. Potential Difference")
#lines(cyanide$t, cyanide$dPsi, cex = 0.01)

pdf("./dataVis/potentialDifferenceCyanide.pdf")

plot(cyanide$t[cyanide$t < 3], cyanide$dPsi[cyanide$t < 3], cex = 0.01,
     xlab = "Time (s)", ylab = "Mem. Potential Difference (mV)", col = "red")
lines(cyanide$t[cyanide$t < 3], cyanide$dPsi[cyanide$t < 3], cex = 0.01,
      col = "red")
points(baselines$t[baselines$t < 3], baselines$dPsi[baselines$t < 3],
       cex = 0.01)
lines(baselines$t[baselines$t < 3], baselines$dPsi[baselines$t < 3],
      cex = 0.01)

dev.off()

pdf("./dataVis/potentialDifferenceCyanideProp.pdf")

plot(cyanide$t[cyanide$t < 3], 
     cyanide$dPsi[cyanide$t < 3]/192.2807, 
     cex = 0.01,
     xlab = "Time (s)", ylab = "Normalized Differnce")
lines(cyanide$t[cyanide$t < 3], 
      cyanide$dPsi[cyanide$t < 3]/192.2807, 
      cex = 0.01)

dev.off()

baselines <- data.table::fread("./results/resultsATP.csv")
cyanide1 <- data.table::fread("./results/resultsCyanide.csv")
cyanide2 <- data.table::fread("./results/resultsCyanide2.csv")
cyanide3 <- data.table::fread("./results/resultsCyanide3.csv")
cyanide10 <- data.table::fread("./results/resultsCyanide10.csv")
hleak10 <- data.table::fread("./results/resultsHLeak9.csv")

pdf("./dataVis/cyanideBinding.pdf")

plot(baselines$t[baselines$t < 45], baselines$dPsi[baselines$t < 45], 
     cex = 0.1, ylim = c(100, 270), ylab = "Potential Gradient", xlab = "Time")
lines(baselines$t[baselines$t < 45], baselines$dPsi[baselines$t < 45], cex = 0.1)
points(cyanide1$t[cyanide1$t < 45], cyanide1$dPsi[cyanide1$t < 45], cex = 0.1, 
       col = "orange")
lines(cyanide1$t[cyanide1$t < 45], cyanide1$dPsi[cyanide1$t < 45], cex = 0.1, 
      col = "orange")
points(cyanide2$t[cyanide2$t < 45], cyanide2$dPsi[cyanide2$t < 45], cex = 0.1, 
       col = "red")
lines(cyanide2$t[cyanide2$t < 45], cyanide2$dPsi[cyanide2$t < 45], cex = 0.1,
      col = "red")
points(cyanide3$t[cyanide3$t < 45], cyanide3$dPsi[cyanide3$t < 45], cex = 0.1, 
       col = "blue")
lines(cyanide3$t[cyanide3$t < 45], cyanide3$dPsi[cyanide3$t < 45], cex = 0.1,
      col = "blue")
points(cyanide10$t[cyanide10$t < 45], cyanide10$dPsi[cyanide10$t < 45], cex = 0.1, 
       col = "cyan")
lines(cyanide10$t[cyanide10$t < 45], cyanide10$dPsi[cyanide10$t < 45], cex = 0.1,
      col = "cyan")
points(hleak10$t[hleak10$t < 45], hleak10$dPsi[hleak10$t < 45], cex = 0.1, 
       col = "violet")
lines(hleak10$t[hleak10$t < 45], hleak10$dPsi[hleak10$t < 45], cex = 0.1,
      col = "violet")

dev.off()

hleak1 <- data.table::fread("./results/resultsHLeak1.csv")
hleak2 <- data.table::fread("./results/resultsHLeak2.csv")
hleak3 <- data.table::fread("./results/resultsHLeak3.csv")
hleak4 <- data.table::fread("./results/resultsHLeak4.csv")
hleak5 <- data.table::fread("./results/resultsHLeak5.csv")
hleak6 <- data.table::fread("./results/resultsHLeak6.csv")
hleak7 <- data.table::fread("./results/resultsHLeak7.csv")
hleak8 <- data.table::fread("./results/resultsHLeak8.csv")

pdf('./dataVis/hleakNormoxic.pdf')
plot(baselines$t[baselines$t < 45 & baselines$t > 10], 
     baselines$dPsi[baselines$t < 45 & baselines$t > 10], 
     cex = 0.1, ylim = c(165, 168), ylab = "Potential Gradient", xlab = "Time")
lines(baselines$t[baselines$t < 45 & baselines$t > 10], 
      baselines$dPsi[baselines$t < 45 & baselines$t > 10], cex = 0.1)
lines(hleak1$t[hleak1$t < 45 & hleak1$t > 10], 
      hleak1$dPsi[hleak1$t < 45 & hleak1$t > 10], cex = 0.1, 
     col = "violet")
lines(hleak2$t[hleak2$t < 45 & hleak2$t > 10], 
      hleak2$dPsi[hleak2$t < 45 & hleak2$t > 10], cex = 0.1, 
      col = "coral")
lines(hleak3$t[hleak3$t < 45 & hleak3$t > 10], 
      hleak3$dPsi[hleak3$t < 45 & hleak3$t > 10], cex = 0.1, 
      col = "blue")
lines(hleak4$t[hleak4$t < 45 & hleak4$t > 10], 
      hleak4$dPsi[hleak4$t < 45 & hleak4$t > 10], cex = 0.1, 
      col = "green")
lines(hleak5$t[hleak5$t < 45 & hleak5$t > 10], 
      hleak5$dPsi[hleak5$t < 45 & hleak5$t > 10], cex = 0.1, 
      col = "springgreen")
lines(hleak6$t[hleak6$t < 45 & hleak6$t > 10], 
      hleak6$dPsi[hleak6$t < 45 & hleak6$t > 10], cex = 0.1, 
      col = "gold")
lines(hleak7$t[hleak7$t < 45 & hleak7$t > 10], 
      hleak7$dPsi[hleak7$t < 45 & hleak7$t > 10], cex = 0.1, 
      col = "orange")
lines(hleak8$t[hleak8$t < 45 & hleak8$t > 10], 
      hleak8$dPsi[hleak8$t < 45 & hleak8$t > 10], cex = 0.1, 
      col = "tomato")
lines(hleak10$t[hleak10$t < 45 & hleak10$t > 10], 
      hleak10$dPsi[hleak10$t < 45 & hleak10$t > 10], cex = 0.1, 
      col = "red")
dev.off()
