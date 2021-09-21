if (!grepl("liverScripts/modelScripts", getwd())) {
  setwd("./liverScripts/modelScripts")
}

hypoxia <- list()
hypoxia[[1]] <- data.table::fread("../results/resultsHypoxiaExtreme.csv")
hypoxia[[1]] <- hypoxia[[1]][hypoxia[[1]]$t > 5000 &
                                   hypoxia[[1]]$t < 25000,]

for (i in 0:9) {
  hypoxia[[i+2]] <- data.table::fread(paste0("../results/resultsHypoxia", 
                                           i, ".csv"))
  hypoxia[[i+2]] <- hypoxia[[i+2]][hypoxia[[i+2]]$t > 5000 &
                                     hypoxia[[i+2]]$t < 25000,]
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
