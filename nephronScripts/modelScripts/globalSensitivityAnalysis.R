if (!grepl("modelScripts", getwd())) {
  setwd("./modelScripts")
}

totalSi <- dir("../results")[grep("totalSi", dir("../results"))]
firstSi <- dir("../results")[grep("firstSi", dir("../results"))]

totalSiList <- list()
firstSiList <- list()

for (i in seq_along(totalSi)) {
  totalSiList[[i]] <- read.csv(paste0("../results/", totalSi[i]))$ST
  firstSiList[[i]] <- read.csv(paste0("../results/", firstSi[i]))$S1
}

paramNames <- read.csv("../results/firstSi0.csv")$X

totalSiList <- do.call(cbind, totalSiList)
firstSiList <- do.call(cbind, firstSiList)

image(abs(firstSiList), zlim = c(0, 1.5))
image(abs(totalSiList), zlim = c(0, 1.5))

importantFirstPars <- paramNames[rowSums(abs(firstSiList)) != 0]
importantTotalPars <- paramNames[rowSums(abs(totalSiList)) != 0]

vimportantFirstPars <- paramNames[rowSums(abs(firstSiList)) > 1]
vimportantTotalPars <- paramNames[rowSums(abs(totalSiList)) > 1]

stateNames <- c("Matrix H", "Electrical Potential", 
                "Matrix ATP", "Matrix ADP",
                "Matrix GTP", "Matrix GDP", "Matrix Pi", 
                "Matrix NADH", "Matrix QH2", "Matrix OAA", 
                "Matrix ACCOA", "Matrix CIT", "Matrix ICIT",
                "Matrix AKG", "Matrix SCOA", "Matrix COASH", 
                "Matrix SUC", "Matrix FUM",
                "Matrix MAL", "Matrix GLU", "Matrix ASP", 
                "Matrix K", "Matrix Mg",
                "Reduced CytC", "IM ATP", 
                "IM ADP", "IM AMP", "IM Pi",
                "Cytosol ATP", "Cytosol ADP", "Cytosol Pi", "Matrix PYR",
                "IM PYR", "IM CIT", "Cytosol CIT", 
                "IM AKG", "Cytosol AKG", "IM SUC", "Cytosol SUC", 
                "IM MAL", "Cytosol MAL", "IM ASP", "Cytosol ASP", 
                "IM GLU", "Cytosol GLU",
                "Cytosol PCr", "Cytosol AMP")
paramNames <- c("PDH Activity", "ICITS Activity", "ACON Activity", 
                "ISOD Activity", "AKGD Activity", "SCOAS Activity", 
                "SDH Activity", "FUM Activity", "MDH Activity", 
                "NDK Activity", "GOT Activity", 
                "PYRH Activity", "GLUH Activity", "CITMAL Activity", 
                "AKGMAL Activity", "MALPI Activity", "ASPGLU Activity", 
                "SUCMAL Activity", "HK Activity", "Complex I Activity", 
                "Complex III Activity", 
                "Complex IV Activity", "ATPase Activity", "ANT Activity", 
                "H/Pi Antiport Activity", "k_PiH", 
                "KH Antiport Activity", "H Leak", "k_Pi1", 
                "k_Pi2", "Mito Water Fraction", "Total CytC", 
                "Total CoQ", "Total NAD/NADH", 
                "Total FAD/FADH2", "k_O2", "k_mADP", 
                "K Leak", "Max. ATP consumption", "Glycolytic ATP")

# stateNames[colSums(abs(firstSiList) > 0.01) > 5]
# stateNames[colSums(abs(totalSiList) > 0.01) > 5]

image(abs(firstSiList[rowSums(abs(firstSiList)) != 0, ]), zlim = c(0, 1.5))
image(abs(totalSiList[rowSums(abs(totalSiList)) != 0, ]), zlim = c(0, 1.5))

neg <- max(abs(firstSiList[firstSiList < 0]))
print(neg)
firstSiList[firstSiList < neg+1e-3] <- 0
firstSiList[firstSiList > 1] <- 1

library(plot.matrix)

colnames(firstSiList) <- colnames(totalSiList) <- stateNames
rownames(firstSiList) <- rownames(totalSiList) <- paramNames

pdf("../dataVis/firstGlobalSensitivity.pdf")
par(mar = c(5.1, 6.1, 2.1, 4.1), las = 2, cex.axis = 0.5)
plot(t(firstSiList), col = viridis::viridis(4), xlab = "", ylab = "", main = "")
dev.off()

pdf("../dataVis/totalGlobalSensitivity.pdf")
par(mar = c(6.1, 6.1, 2.1, 4.1), las = 2, cex.axis = 0.5)
plot(t(totalSiList), col = viridis::viridis, xlab = "", ylab = "", main = "")
dev.off()

## Important states
pdf("../dataVis/firstGlobalSensitivityImportant.pdf")
par(mar = c(10.1, 8.1, 2.1, 4.1), las = 2, cex.axis = 1)
plot(t(firstSiList[rowSums(abs(firstSiList[ , c(2, 8, 9, 24, 29)])) > 0.25
  , c(2, 8, 9, 24, 29)]), col = viridis::viridis(4), xlab = "", ylab = "", main = "")
dev.off()

pdf("../dataVis/totalGlobalSensitivityImportant.pdf")
par(mar = c(10.1, 8.1, 2.1, 4.1), las = 2, cex.axis = 1)
plot(t(totalSiList[rowSums(totalSiList[ , c(2, 8, 9, 24, 29)]) > 0.25
                   , c(2, 8, 9, 24, 29)]), col = viridis::viridis, xlab = "", 
     ylab = "", main = "")
dev.off()
