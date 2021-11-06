mTALfit <- list()
default <- data.table::fread("./results/resultsmTAL.csv")

for (i in 0:6) {
  mTALfit[[i+1]] <- data.table::fread(paste0("./results/resultsmTALfit", i, ".csv"))
}

plot(x = NULL, ylim = c(0, 0.0025), xlim = c(1, 3e1), ylab = "ATP", xlab = "Time")
for (i in 1:7) {
  lines(mTALfit[[i]]$t, mTALfit[[i]]$ATP_c, col = c("black", "red")[(i == 5)+1])
}
lines(default$t, default$ATP_c, col = "blue")

plot(x = NULL, ylim = c(140, 170), xlim = c(1, 3e1), ylab = "dPsi", xlab = "Time")
for (i in 1:7) {
  lines(mTALfit[[i]]$t, mTALfit[[i]]$dPsi, col = c("black", "red")[(i == 5)+1])
}
lines(default$t, default$dPsi, col = "blue")

plot(x = NULL, ylim = c(0, 1), xlim = c(1, 3e1), ylab = "Redox state", xlab = "Time")
nadtotPerm <- 0.824e-4
for (i in 1:7) {
  nadtotTemp <- nadtotPerm
  if (i == 7) {
    nadtotTemp <- 10*nadtotTemp
  }
  lines(mTALfit[[i]]$t, mTALfit[[i]]$NADH_x/(nadtotTemp), 
        col = c("black", "red")[(i == 5)+1])
}
lines(default$t, default$NADH_x/(nadtotPerm), col = "blue")


## Mitochondrial water fraction appears to not matter at all (one less param to fit)
## NADTOT does appear to change things though - which could be good because I was
## hoping to fit redox state
## Increasing ylim will reveal that low NADTOT because of the extremely imbalanced
## redox state can push the dPsi really high in the transient