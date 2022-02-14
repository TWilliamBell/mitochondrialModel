if (!grepl("modelScripts", getwd())) {
  setwd("modelScripts")
}

## In order to run this file, you should first run localSensitivity.py (and if you want to make the
## combined plot, then you should also run localSensitivitymTAL.py).

localSensitivity <- read.csv("../results/localSensitivityStates.csv")
localSensitivity$X <- NULL
dirLS <- matrix(0, nrow = 40, ncol = 63)

typical <- c(4.34260851e-08, 1.55176612e+02, 9.84096832e-04, 1.39090317e-03,
             3.75000000e-04, 2.14730531e-03, 2.85272469e-03, 7.78911327e-03,
             6.15169703e-05, 2.84209686e-04, 1.49368365e-10, 3.00851364e-03,
             8.63430185e-04, 1.70942773e-05, 1.18871334e-03, 6.58407994e-08,
             1.28131386e-06, 2.04542575e-05, 2.28931607e-05, 1.08154743e-05,
             8.54706553e-03, 7.46929765e-10, 1.29102218e-01, 1.95937766e-03,
             6.69000000e-05, 2.14000000e-02, 6.48462171e-04, 2.27659478e-03,
             4.20507320e-04, 5.10751058e-05, 3.73659863e-03, 6.30960000e-08,
             1.00300000e-03, 1.30000000e-01, 2.26959646e-03, 4.27505639e-04,
             3.73841777e-03, 6.30960000e-08, 4.00000000e-04, 1.10000000e-01,
             1.99840050e-04, 1.53093432e-04, 1.53750000e-04, 3.91272198e-04,
             3.91272203e-04, 5.35814089e-04, 5.35814086e-04, 9.32724461e-06,
             9.32724461e-06, 4.87510034e-06, 4.87510034e-06, 1.90908103e-07,
             1.90908104e-07, 5.34736887e-03, 5.34736887e-03, 1.50000000e-05,
             6.25000000e-05, 7.08405956e-05, 7.08405956e-05, 7.25000000e-03,
             9.12500000e-05, 4.08102233e-03, 5.29520697e-05)

for (i in 0:38) {
  dirLS[i+1, ] <- 100*(unlist(localSensitivity[2*i+1, ])-
                     unlist(localSensitivity[2*i+2, ]))/typical
}

library('plot.matrix')
image(abs(dirLS))
colnames(dirLS) <- c("H_x", "dPsi", "ATP_x", "ADP_x",
                     "AMP_x", "GTP_x", "GDP_x", "Pi_x", "NADH_x",
                     "QH2_x", "OAA_x", "ACCOA_x", "CIT_x", "ICIT_x",
                     "AKG_x", "SCOA_x", "COASH_x", "SUC_x", "FUM_x",
                     "MAL_x", "GLU_x", "ASP_x", "K_x", "Mg_x",
                     "O2_x", "CO2tot_x", "Cred_i", "ATP_i", "ADP_i",
                     "AMP_i", "Pi_i", "H_i", "Mg_i", "K_i", "ATP_c",
                     "ADP_c", "Pi_c", "H_c", "Mg_c", "K_c", "PYR_x",
                     "PYR_i", "PYR_c", "CIT_i", "CIT_c", "AKG_i",
                     "AKG_c", "SUC_i", "SUC_c", "MAL_i", "MAL_c",
                     "ASP_i", "ASP_c", "GLU_i", "GLU_c", "FUM_i",
                     "FUM_c", "ICIT_i", "ICIT_c", "GLC_c", "G6P_c",
                     "PCr_c", "AMP_c")
rownames(dirLS) <- c("Vmax1", "Vmax2", "Vmax3", "Vmax4", "Vmax5", "Vmax6", "Vmax7", 
                     "Vmax8", "Vmax9", "Vmax10", "Vmax11", "x_PYR_H", "x_GLU_H", "x_CIT_MAL", 
                     "x_AKG_MAL", "x_MAL_PI", "x_ASP_GLU", "x_SUC_MAL", "X_HK", "x_C1", "x_C3", 
                     "x_C4", "x_F1", "x_ANT", "x_Pi1", "k_PiH", "x_KH", "x_Hle", "k_Pi1", 
                     "k_Pi2", "W_m", "Ctot", "Qtot", "NADtot", "FADtot", "k_O2", "k_mADP", 
                     "pW", "J_AtC", "glyc")
dirLS <- t(dirLS)

## Figure 3.1
pdf("../dataVis/localSensitivityOverall.pdf", width = 20, height = 20)
par(mar = c(10.1, 10.1, 2.1, 4.1), las = 2)
tempDirLS <- dirLS
tempDirLS <- t(tempDirLS)
rownames(tempDirLS) <- c("PDH Activity", "ICITS Activity", "ACON Activity", 
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
colnames(tempDirLS) <- c("Matrix H", "Electrical Potential", 
                         "Matrix ATP", "Matrix ADP",
                         "Matrix AMP", "Matrix GTP", "Matrix GDP", "Matrix Pi", 
                         "Matrix NADH", "Matrix QH2", "Matrix OAA", 
                         "Matrix ACCOA", "Matrix CIT", "Matrix ICIT",
                         "Matrix AKG", "Matrix SCOA", "Matrix COASH", 
                         "Matrix SUC", "Matrix FUM",
                         "Matrix MAL", "Matrix GLU", "Matrix ASP", 
                         "Matrix K", "Matrix Mg",
                         "Oxygen", "Carbon Dioxide", "Reduced CytC", "IM ATP", 
                         "IM ADP", "IM AMP", "IM Pi", "IM H", "IM Mg", "IM K", 
                         "Cytosol ATP", "Cytosol ADP", "Cytosol Pi", 
                         "Cytosol H", "Cytosol Mg", "Cytosol K", "Matrix PYR",
                         "IM PYR", "Cytosol PYR", "IM CIT", "Cytosol CIT", 
                         "IM AKG", "Cytosol AKG", "IM SUC", "Cytosol SUC", 
                         "IM MAL", "Cytosol MAL", "IM ASP", "Cytosol ASP", 
                         "IM GLU", "Cytosol GLU", "IM FUM",
                         "Cytosol FUM", "IM ICIT", "Cytosol ICIT", 
                         "Cytosol GLC", "Cytosol G6P",
                         "Cytosol PCr", "Cytosol AMP")
tempDirLS <- t(tempDirLS)
plot(tempDirLS, col = viridis::viridis, xlab = "", ylab = "", main = "")
dev.off()

importantStatedirLS <- dirLS[rowSums(abs(dirLS)) != 0, ]
image(abs(importantStatedirLS))
plot(abs(importantStatedirLS), col = viridis::viridis)

importantParamdirLS <- dirLS[ , colSums(abs(dirLS)) != 0]
image(abs(importantParamdirLS))
plot(importantParamdirLS, col = viridis::viridis)

vimportantStatedirLS <- dirLS[rowSums(abs(dirLS)) > 1, ]
image(abs(importantStatedirLS))
plot(vimportantStatedirLS, col = viridis::viridis)

vimportantParamdirLS <- dirLS[ , colSums(abs(dirLS)) > 1]
image(abs(importantParamdirLS))
plot(vimportantParamdirLS, col = viridis::viridis)

bothImportantdirLS <- dirLS[rowSums(abs(dirLS)) != 0, colSums(abs(dirLS)) != 0]
image(abs(bothImportantdirLS))
plot(bothImportantdirLS, col = viridis::viridis)

bothvimportantdirLS <- dirLS[rowSums(abs(dirLS)) > 1, colSums(abs(dirLS)) > 1]
image(abs(bothvimportantdirLS))
par(mar = c(4.1, 4.1, 2.1, 3.1), cex.axis = 0.5, las = 2)
plot(bothvimportantdirLS, col = viridis::viridis, xlab = "", ylab = "", main = "")


paramNames <- c("Vmax1", "Vmax2", "Vmax3", "Vmax4", "Vmax5", "Vmax6", "Vmax7", 
                "Vmax8", "Vmax9", "Vmax10", "Vmax11", "x_PYR_H", "x_GLU_H", "x_CIT_MAL", 
                "x_AKG_MAL", "x_MAL_PI", "x_ASP_GLU", "x_SUC_MAL", "X_HK", "x_C1", "x_C3", 
                "x_C4", "x_F1", "x_ANT", "x_Pi1", "k_PiH", "x_KH", "x_Hle", "k_Pi1", 
                "k_Pi2", "W_m", "Ctot", "Qtot", "NADtot", "FADtot", "k_O2", "k_mADP", 
                "pW", "J_AtC", "glyc")
importantParams <- paramNames[colSums(abs(dirLS)) > 1]

stateNames <- c("H_x", "dPsi", "ATP_x", "ADP_x",
                  "AMP_x", "GTP_x", "GDP_x", "Pi_x", "NADH_x",
                  "QH2_x", "OAA_x", "ACCOA_x", "CIT_x", "ICIT_x",
                  "AKG_x", "SCOA_x", "COASH_x", "SUC_x", "FUM_x",
                  "MAL_x", "GLU_x", "ASP_x", "K_x", "Mg_x",
                  "O2_x", "CO2tot_x", "Cred_i", "ATP_i", "ADP_i",
                  "AMP_i", "Pi_i", "H_i", "Mg_i", "K_i", "ATP_c",
                  "ADP_c", "Pi_c", "H_c", "Mg_c", "K_c", "PYR_x",
                  "PYR_i", "PYR_c", "CIT_i", "CIT_c", "AKG_i",
                  "AKG_c", "SUC_i", "SUC_c", "MAL_i", "MAL_c",
                  "ASP_i", "ASP_c", "GLU_i", "GLU_c", "FUM_i",
                  "FUM_c", "ICIT_i", "ICIT_c", "GLC_c", "G6P_c",
                  "PCr_c", "AMP_c")
importantStates <- stateNames[rowSums(abs(dirLS)) > 1]

par(mar = c(5.1, 4.1, 4.1, 2.1), cex.axis = 0.5)
barplot(dirLS[35, ]) ## ATP_c
barplot(dirLS[2, ]) ## dPsi
barplot(dirLS[9, ]) ## NADH_x
barplot(dirLS[10, ]) ## QH2_x
barplot(dirLS[27, ]) ## Cred_i

dirLSsmall <- dirLS[c(2, 9, 10, 27, 35), colSums(abs(dirLS[ 
                         c(2, 9, 10, 27, 35), ])) > 0.5]
colnames(dirLSsmall) <- c("PDH Activity", "AKGD Activity", 
                          "Complex III Activity", "Complex IV Activity",
                          "ANT Activity", "k_Pi1", "k_Pi2",
                          "Total CytC", "Total CoQ", "Total NAD+/NADH", 
                          "k_O2", "Max. ATP Consumption")
rownames(dirLSsmall) <- c("Electrical Potential", "NADH", "QH2", 
                          "Reduced CytC", "Cytosolic ATP")

pdf("../dataVis/importantLocalSensitivity.pdf", width = 10)
par(mar = c(10.1, 8.1, 1.1, 4.1), las = 2, cex.axis = 1)
plot(dirLSsmall, 
     col = viridis::viridis(12), main = "",
     xlab = "", ylab = "")
dev.off()

dirLSPT <- dirLS
## Combined figure
# To run this code: 1. Source localSensitivityAnalysis.R,
# 2. Source localSensitivityAnalysismTAL.R,
# 3. Then run the commented code
# 
# dirLSsmallPT <- dirLSPT[c(2, 9, 10, 27, 35), colSums(abs(dirLSPT[
#   c(2, 9, 10, 27, 35), ])) > 0.5 | colSums(abs(dirLS[
#     c(2, 9, 10, 27, 35), ])) > 0.5]
# dirLSsmallTAL <- dirLS[c(2, 9, 10, 27, 35), colSums(abs(dirLSPT[
#   c(2, 9, 10, 27, 35), ])) > 0.5 | colSums(abs(dirLS[
#     c(2, 9, 10, 27, 35), ])) > 0.5]
# 
# colnames(dirLSsmallTAL) <- colnames(dirLSsmallPT) <- c("PDH Activity", "AKGD Activity",
#                           "Complex II Activity",
#                           "Complex III Activity", "Complex IV Activity",
#                           "ANT Activity", "k_Pi1", "k_Pi2",
#                           "Total CytC", "Total CoQ", "Total NAD+/NADH",
#                           "Complex IV Oxygen Affinity",
#                           "Max. ATP Consumption")
# rownames(dirLSsmallTAL) <- rownames(dirLSsmallPT) <- c("Electrical Potential", "NADH", "QH2",
#                           "Reduced CytC", "Cytosolic ATP")
# 
# ## Figure 3.3
# 
# pdf("../dataVis/combinedLS.pdf", width = 10)
# par(mar = c(12.1, 8.1, 1.1, 2.1), las = 2, cex.axis = 1, mfrow = c(1, 2))
# plot(dirLSsmallPT,
#      col = viridis::viridis, main = "",
#      xlab = "", ylab = "", breaks = seq(from = -1.5, to = 1.5, by = 0.25),
#      key = NULL)
# par(mar = c(12.1, 6.1, 1.1, 4.1))
# plot(dirLSsmallTAL,
#      col = viridis::viridis, main = "",
#      xlab = "", ylab = "", breaks = seq(from = -1.5, to = 1.5, by = 0.25))
# dev.off()
# 
