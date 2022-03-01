## Before running this script, run hypoxia.py

if (!grepl("liverScripts/modelScripts", getwd())) {
  setwd("./liverScripts/modelScripts")
}

hypoxia <- list()
hypoxia[[1]] <- data.table::fread("../results/resultsHypoxiaExtreme.csv")
hypoxia[[1]] <- hypoxia[[1]]

for (i in 0:9) {
  hypoxia[[i+2]] <- data.table::fread(paste0("../results/resultsHypoxia", 
                                           i, ".csv"))
}

colVal <- rainbow(11)
pdf("../dataVis/hypoxiaVis.pdf")
plot(NULL,#hypoxia[[1]]$t, hypoxia[[1]]$ATP_c, cex = 0.001, col = colVal[1], 
     xlim = c(5000, 25000),
     ylim = c(0, 7), xlab = "Time (s)", ylab = "Cytosolic ATP (mM)")

for (i in 1:11) {
  lines(hypoxia[[i]]$t, hypoxia[[i]]$ATP_c*1000, col = colVal[i])
}
dev.off()

atpLevels <- sapply(hypoxia, function(x) min(x$ATP_c))*1000
o2Levels <- c(0.5/36, (1:10)/10)

## Figure 4.4
pdf("../dataVis/hypoxiaResponseLiver.pdf")
par(cex.axis = 2, cex.lab = 1.5)
plot(o2Levels, atpLevels, col = "red", cex = 0.5, xlab = "Fold-Change in Oxygen Tension",
     ylab = "Cytosolic ATP (mM)")
lines(o2Levels, atpLevels, col = "red")
dev.off()
