## Capacitance

cyanKBif <- list()
cyanKHBif <- list()
cyanKHHBif <- list()

cyanKBif[[1]] <- cyanKHBif[[1]] <- data.table::fread("./results/resultsCyanide2.csv")
cyanKHHBif[[1]] <- data.table::fread("./results/resultsHLeak9.csv")

for (i in 2:10) {
  cyanKBif[[i]] <- data.table::fread(paste0("./results/resultsKCyan", i, ".csv"))
  cyanKHBif[[i]] <- data.table::fread(
                      paste0("./results/resultsKHAntiporterCyan", i, ".csv"))
  cyanKHHBif[[i]] <- data.table::fread(
                      paste0("./results/resultsKHAntiporterHleakCyan", i, ".csv"))
}

plot(cyanKBif[[1]]$t[cyanKBif[[1]]$t < 300], 
     cyanKBif[[1]]$dPsi[cyanKBif[[1]]$t < 300], cex = 0.1, col = "red", 
     ylim = c(80, 150), ylab = "dPsi", xlab = "t")
lines(cyanKBif[[1]]$t[cyanKBif[[1]]$t < 300], 
  cyanKBif[[1]]$dPsi[cyanKBif[[1]]$t < 300], cex = 0.1, col = "red")
lines(cyanKHBif[[1]]$t[cyanKHBif[[1]]$t < 300], 
      cyanKHBif[[1]]$dPsi[cyanKHBif[[1]]$t < 300], cex = 0.1, col = "blue")
lines(cyanKBif[[9]]$t[cyanKBif[[9]]$t < 300], 
      cyanKBif[[9]]$dPsi[cyanKBif[[9]]$t < 300], cex = 0.1, col = "green")
lines(cyanKHBif[[9]]$t[cyanKHBif[[9]]$t < 300], 
      cyanKHBif[[9]]$dPsi[cyanKHBif[[9]]$t < 300], cex = 0.1, col = "orange")
lines(cyanKHHBif[[9]]$t[cyanKHHBif[[9]]$t < 300], 
      cyanKHHBif[[9]]$dPsi[cyanKHHBif[[9]]$t < 300], cex = 0.1, col = "springgreen")

## Maxxed potassium leak, hydrogen concentrations, cyanide poisoning
tail(cyanKBif[[10]]$H_x, 1)
tail(cyanKBif[[10]]$H_i, 1)

## Reduced K-H exchange, hydrogen concentrations, cyanide poisoning
tail(cyanKHBif[[10]]$H_x, 1)
tail(cyanKHBif[[10]]$H_i, 1)

## Just cyanide, hydrogen concentrations
tail(cyanKBif[[1]]$H_x, 1)
tail(cyanKBif[[1]]$H_i, 1)
