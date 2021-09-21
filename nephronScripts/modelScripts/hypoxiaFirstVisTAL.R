hypoxia <- list()

for (i in 0:9) {
  hypoxia[[i+1]] <- data.table::fread(
    paste0("../results/resultsHypoxiaTAL", i, ".csv"))
}

a <- rainbow(10)

pdf("../dataVis/hypoxiaTALparam.pdf")
par(cex.lab = 1.5, cex.axis = 1.5)
for (i in 1:10) {
  if (i == 1) {
    plot(hypoxia[[i]]$t[hypoxia[[i]]$t > 5000 & hypoxia[[i]]$t < 25000], 
         hypoxia[[i]]$ATP_c[hypoxia[[i]]$t > 5000 & hypoxia[[i]]$t < 25000], 
         ylab = "[ATP]_c (M)",
         xlab = "Time (s)", cex = 0.1, col = a[i], 
         ylim = c(0.00, 0.0030))
    lines(hypoxia[[i]]$t[hypoxia[[i]]$t > 5000 & hypoxia[[i]]$t < 25000],
          hypoxia[[i]]$ATP_c[hypoxia[[i]]$t > 5000 & hypoxia[[i]]$t < 25000],
          ylab = "[ATP]_c",
          xlab = "Time", cex = 0.1, col = a[i],
          ylim = c(0.00, 0.0030))
  } else {
    points(hypoxia[[i]]$t[hypoxia[[i]]$t > 5000 & hypoxia[[i]]$t < 25000],
           hypoxia[[i]]$ATP_c[hypoxia[[i]]$t > 5000 & hypoxia[[i]]$t < 25000],
           ylab = "[ATP]_c",
           xlab = "Time", cex = 0.1, col = a[i],
           ylim = c(0.00, 0.0030))
    lines(hypoxia[[i]]$t[hypoxia[[i]]$t > 5000 & hypoxia[[i]]$t < 25000],
          hypoxia[[i]]$ATP_c[hypoxia[[i]]$t > 5000 & hypoxia[[i]]$t < 25000],
          ylab = "[ATP]_c",
          xlab = "Time", cex = 0.1, col = a[i],
          ylim = c(0.00, 0.0030))
  }
}

b <- as.character(((1:10)/10)*10)
b[1] <- paste(b[1], "mmHg")

legend(x = c(20800, 25300), y = c(0.0008, 0.0022), legend = b, fill = a,
       cex = 1)

dev.off()
