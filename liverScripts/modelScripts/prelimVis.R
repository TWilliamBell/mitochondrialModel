#setwd("C:/Users/twill/Desktop/Folder/School/MMath Thesis/mitochondrialModel/liverScripts/modelScripts")

results <- data.table::fread("../results/resultsATP.csv")

plot(results$t, results$dPsi, cex = 0.1, col = "red", ylab = "dPsi", xlab = "t")
lines(results$t, results$dPsi, col = "red")

plot(results$t, results$ATP_c, cex = 0.1, col = "red", ylab = "ATP_c", xlab = "t")
lines(results$t, results$ATP_c, col = "red")
