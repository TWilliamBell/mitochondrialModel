nigericin <- data.table::fread("./results/resultsNigericin.csv")

plot(nigericin$t, nigericin$dPsi, cex = 0.1, ylim = c(167, 168))
lines(nigericin$t, nigericin$dPsi)

norm <- data.table::fread("./results/resultsATP.csv")

plot(norm$t, norm$dPsi, cex = 0.1, ylim = c(167, 168))
lines(norm$t, norm$dPsi)

print(tail(nigericin$dPsi, 1)/tail(norm$dPsi, 1))
