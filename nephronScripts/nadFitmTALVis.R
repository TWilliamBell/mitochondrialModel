mTALfit <- list()
default <- data.table::fread("./results/resultsmTAL.csv")

for (i in 0:10) {
  mTALfit[[i+1]] <- data.table::fread(paste0("./results/resultsmTALNAD", i, ".csv"))
}

plot(x = NULL, ylim = c(0, 0.0025), xlim = c(1, 3e1), ylab = "ATP", xlab = "Time")
for (i in 1:10) {
  lines(mTALfit[[i]]$t, mTALfit[[i]]$ATP_c, col = c("black", "red")[(i == 10)+1])
}
lines(default$t, default$ATP_c, col = "blue")

plot(x = NULL, ylim = c(140, 170), xlim = c(1, 3e1), ylab = "dPsi", xlab = "Time")
for (i in 1:10) {
  lines(mTALfit[[i]]$t, mTALfit[[i]]$dPsi, col = c("black", "red")[(i == 10)+1])
}
lines(default$t, default$dPsi, col = "blue")

plot(x = NULL, ylim = c(0, 1), xlim = c(1, 3e1), ylab = "Redox state", xlab = "Time")
nadtotPerm <- 0.824e-3
for (i in 1:10) {
  nadtotTemp <- nadtotPerm/(19+i)
  lines(mTALfit[[i]]$t, mTALfit[[i]]$NADH_x/(nadtotTemp), 
        col = c("black", "red")[(i == 10)+1])
}
lines(default$t, default$NADH_x/(nadtotPerm), col = "blue")
