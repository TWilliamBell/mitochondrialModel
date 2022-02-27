## In order to run this file, first run hypoxiaTAL.py

hypoxia <- list()

for (i in 0:9) {
  hypoxia[[i+1]] <- data.table::fread(
    paste0("../results/resultsHypoxiaTAL", i, ".csv"))
  hypoxia[[i+1]]$t <- (hypoxia[[i+1]]$t/3600)-5000/3600
}

a <- rainbow(10)

## Figure 3.19
pdf("../dataVis/hypoxiaTALparam.pdf")
par(cex.axis = 1.5, cex.lab = 1.5)
for (i in 1:10) {
  if (i == 1) {
    plot(hypoxia[[i]]$t[hypoxia[[i]]$t > 0 & hypoxia[[i]]$t < 20000/3600], 
         hypoxia[[i]]$ATP_c[hypoxia[[i]]$t > 0 & hypoxia[[i]]$t < 20000/3600]*1000, 
         ylab = "Cytosolic ATP (mM)",
         xlab = "Time (hrs)", cex = 0.1, col = a[i], 
         ylim = c(0.00, 0.0026)*1000)
    lines(hypoxia[[i]]$t[hypoxia[[i]]$t > 0 & hypoxia[[i]]$t < 20000/3600],
          hypoxia[[i]]$ATP_c[hypoxia[[i]]$t > 0 & hypoxia[[i]]$t < 20000/3600]*1000,
          ylab = "[ATP]_c",
          xlab = "Time", cex = 0.1, col = a[i], add = T,
          ylim = c(0.00, 0.0026)*1000)
  } else {
    points(hypoxia[[i]]$t[hypoxia[[i]]$t > 0 & hypoxia[[i]]$t < 20000/3600],
           hypoxia[[i]]$ATP_c[hypoxia[[i]]$t > 0 & hypoxia[[i]]$t < 20000/3600]*1000,
           ylab = "Cytosolic ATP (mM)",
           xlab = "Time", cex = 0.1, col = a[i], add = T,
           ylim = c(0.00, 0.0026)*1000)
    lines(hypoxia[[i]]$t[hypoxia[[i]]$t > 0 & hypoxia[[i]]$t < 20000/3600],
          hypoxia[[i]]$ATP_c[hypoxia[[i]]$t > 0 & hypoxia[[i]]$t < 20000/3600]*1000,
          ylab = "[ATP]_c",
          xlab = "Time", cex = 0.1, col = a[i], add = T,
          ylim = c(0.00, 0.0026)*1000)
  }
}

b <- as.character(((1:10)/10)*10)
b[1] <- paste(b[1], "mmHg")

legend(x = c(20800-5000, 25300-5000)/3600, y = c(0.001*1000, 0.0022*1000), legend = b, fill = a,
       cex = 1)

dev.off()


f <- function(x) {
  x <- x[x$t > 1, ]
  min(x$ATP_c)
}

hypoxiaMins <- sapply(hypoxia, f)

o2fn <- splinefun((1:10)/10, 1000*hypoxiaMins)

## Figure 3.20
pdf("../dataVis/hypoxiamTALresponse.pdf")
par(cex.axis = 1.5, cex.lab = 1.5)
plot((1:10)/10, hypoxiaMins*1000, xlab = "Fold Change in Oxygen Tension", ylab = "Cytosolic ATP (mM)")
curve(o2fn, from = 0, to = 1, add = T)
dev.off()
