## Capacitance

normCapBif <- list() ## No cyanide, range of values for capacitance
cyanCapBif <- list() ## Cyanide and a range of values for capacitances

normCapBif[[1]] <- feather::read_feather("./results/resultsATP.feather")
cyanCapBif[[1]] <- data.table::fread("./results/resultsCyanide2.csv")


for (i in 2:9) {
  normCapBif[[i]] <- feather::read_feather(paste0("./results/resultsCap", i, ".feather"))
  cyanCapBif[[i]] <- feather::read_feather(
                      paste0("./results/resultsCapCyan", i, ".feather"))
}

plot(normCapBif[[1]]$t[normCapBif[[1]]$t < 1], 
     normCapBif[[1]]$dPsi[normCapBif[[1]]$t < 1], cex = 0.1, col = "red", 
     ylim = c(100, 350), ylab = "dPsi", xlab = "t")
lines(normCapBif[[1]]$t[normCapBif[[1]]$t < 1], 
  normCapBif[[1]]$dPsi[normCapBif[[1]]$t < 1], cex = 0.1, col = "red")
lines(cyanCapBif[[1]]$t[cyanCapBif[[1]]$t < 1], 
      cyanCapBif[[1]]$dPsi[cyanCapBif[[1]]$t < 1], cex = 0.1, col = "blue")
lines(normCapBif[[9]]$t[normCapBif[[9]]$t < 1], 
      normCapBif[[9]]$dPsi[normCapBif[[9]]$t < 1], cex = 0.1, col = "green")
lines(cyanCapBif[[9]]$t[cyanCapBif[[9]]$t < 1], 
      cyanCapBif[[9]]$dPsi[cyanCapBif[[9]]$t < 1], cex = 0.1, col = "orange")
