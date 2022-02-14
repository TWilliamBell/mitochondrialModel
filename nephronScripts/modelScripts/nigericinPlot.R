if (!grepl("mitochondrialModel/modelScripts", getwd())) {
  setwd("./modelScripts")
}

nigericin <- list()

nigericin[[1]] <- data.table::fread("../results/resultsATP.csv")

nigericin[[2]] <- data.table::fread(paste0("../results/resultsNigericinExp.csv"))


t1 <- 20
t2 <- 60

pdf("./dataVis/nigericin.pdf")
plot(nigericin[[1]]$t[nigericin[[1]]$t > t1 & nigericin[[1]]$t < t2], 
     nigericin[[1]]$dPsi[nigericin[[1]]$t > t1 & nigericin[[1]]$t < t2], 
     ylim = c(0, 200), xlim = c(t1, t2), cex = 0.1, col = "red",
     ylab = "Electric Potential Gradient", xlab = "Time")
lines(nigericin[[1]]$t[nigericin[[1]]$t > t1 & nigericin[[1]]$t < t2], 
      nigericin[[1]]$dPsi[nigericin[[1]]$t > t1 & nigericin[[1]]$t < t2], 
      col = "red", cex = 0.1)
lines(nigericin[[2]]$t[nigericin[[2]]$t > t1 & nigericin[[2]]$t < t2], 
      nigericin[[2]]$dPsi[nigericin[[2]]$t > t1 & nigericin[[2]]$t < t2], 
      cex = 0.01, col = "blue")
arrows(x0 = 40, y0 = 120, y1 = 150, length = 0.15)
text(40, 112, "Nigericin Added")
dev.off()

pdf("./dataVis/nigericinBar.pdf")
barplot(c(159.1, 168.3), names.arg = c("No Nigericin", "Nigericin"),
        ylim = c(0, 200), col = c("turquoise", "tomato"),
        ylab = "Electrical Potential Gradient (mV)")
dev.off()

pdf("./dataVis/nigericinATP.pdf")
plot(nigericin[[1]]$t[nigericin[[1]]$t > t1 & nigericin[[1]]$t < t2], 
     nigericin[[1]]$ATP_c[nigericin[[1]]$t > t1 & nigericin[[1]]$t < t2], 
     ylim = c(2e-3, 3e-3), xlim = c(t1, t2), cex = 0.1, col = "red",
     ylab = "Electric Potential Gradient", xlab = "Time")
lines(nigericin[[1]]$t[nigericin[[1]]$t > t1 & nigericin[[1]]$t < t2], 
      nigericin[[1]]$ATP_c[nigericin[[1]]$t > t1 & nigericin[[1]]$t < t2], 
      col = "red", cex = 0.1)
lines(nigericin[[2]]$t[nigericin[[2]]$t > t1 & nigericin[[2]]$t < t2], 
      nigericin[[2]]$ATP_c[nigericin[[2]]$t > t1 & nigericin[[2]]$t < t2], 
      cex = 0.01, col = "blue")
dev.off()