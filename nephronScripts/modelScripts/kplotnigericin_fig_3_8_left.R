if (!grepl("modelScripts", getwd())) {
  setwd("./modelScripts")
}

## In order to run this file, you first need to run potassium.py

readfn <- function(x) lapply(x, function(x) {data.table::fread(x)})
lowK <- paste0("../results/resultsNigericinKLeak", 1:10, "K1.csv")
normalK <- paste0("../results/resultsNigericin", 1:10, ".csv")
highK <- paste0("../results/resultsNigericinKLeak", 1:10, "K3.csv")

lowKDF <- readfn(lowK)
normalKDF <- readfn(normalK)
highKDF <- readfn(highK)

beforefn <- function(x) sapply(x, function(x) tail(x[x$t < 30, ], 1)$dPsi)

beforelowKtails <- beforefn(lowKDF)
beforenormalKtails <- beforefn(normalKDF)
beforehighKtails <- beforefn(highKDF)

afterfn <- function(x) sapply(x, function(x) tail(x, 1)$dPsi)

afterlowKtails <- afterfn(lowKDF)
afternormalKtails <- afterfn(normalKDF)
afterhighKtails <- afterfn(highKDF)

changelowK <- c(1, afterlowKtails/beforelowKtails)
changenormalK <- c(1, afternormalKtails/beforenormalKtails)
changehighK <- c(1, afterhighKtails/beforehighKtails)

pdf("../dataVis/kleakKHAntiporter.pdf")
par(cex.axis = 1.5, cex.lab = 1.5)
fHigh <- splinefun(x = c(1, 1:10*2), y = changehighK)
curve(fHigh, from = 1, to = 20, col = "red", ylab = "Fold Change in Potential Gradient",
      xlab = "Fold Change in Potassium-Hydrogen Antiporter Activity")
points(x = c(1, 1:10*2), y = changehighK, col = "red")
fNormal <- splinefun(x = c(1, 1:10*2), y = changenormalK)
curve(fNormal, from = 1, to = 20, add = T, col = "green")
points(x = c(1, 1:10*2), y = changenormalK, col = "green")
fLow <- splinefun(x = c(1, 1:10*2), y = changelowK)
curve(fLow, from = 1, to = 20, add = T, col = "blue")
points(x = c(1, 1:10*2), y = changelowK, col = "blue")
dev.off()

## Figure 3.8 first plate
pdf("../dataVis/kleakKHAntiporterAbs.pdf")
par(cex.axis = 1.5, cex.lab = 1.5)
fHigh <- splinefun(x = c(1, 1:10*2), y = c(beforehighKtails[1], afterhighKtails))
curve(fHigh, from = 1, to = 20, col = "red", ylab = "Electrical Potential Gradient (mV)",
      xlab = "Fold Change in Potassium-Hydrogen Antiporter Activity", ylim = c(145, 175))
points(x = c(1, 1:10*2), y = c(beforehighKtails[1], afterhighKtails), col = "red")
fNormal <- splinefun(x = c(1, 1:10*2), y = c(beforenormalKtails[1], afternormalKtails))
curve(fNormal, from = 1, to = 20, add = T, col = "green")
points(x = c(1, 1:10*2), y = c(beforenormalKtails[1], afternormalKtails), col = "green")
fLow <- splinefun(x = c(1, 1:10*2), y = c(beforelowKtails[1], afterlowKtails))
curve(fLow, from = 1, to = 20, add = T, col = "blue")
points(x = c(1, 1:10*2), y = c(beforelowKtails[1], afterlowKtails), col = "blue")
dev.off()
