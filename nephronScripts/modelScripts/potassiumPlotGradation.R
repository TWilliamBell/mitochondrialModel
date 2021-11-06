if (!grepl("modelScripts", getwd())) {
  setwd("modelScripts")
}

potassiumNames <- paste0("../results/resultsNigericin", 1:10, ".feather")

potassium <- lapply(potassiumNames, feather::read_feather)
potassium <- sapply(potassium, function(x) tail(x$dPsi, 1)/159.7)

g <- splinefun(c(1, 1:10*2), c(1, potassium))

pdf("../dataVis/potassiumGradation.pdf", width = 9)
par(cex.lab = 2, cex.axis = 2, mar = c(6.1, 5, 5, 2.1))
plot(c(1, 1:10*2), c(1, potassium), xlab = "Fold Change in K+/H+ Antiport Activity",
     ylab = "Change in Electrical Potential Gradient", cex = 0.5, col = "red", 
     xlim = c(0, 20))
curve(g, from = 0, to = 20, col = "red", add = T, n = 1000)
dev.off()
