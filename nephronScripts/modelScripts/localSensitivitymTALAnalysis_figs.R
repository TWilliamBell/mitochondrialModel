if (!grepl("modelScripts", getwd())) {
  setwd("modelScripts")
}

## In order to run this file, you should first run localSensitivitymTAL.py.

localSensitivity <- read.csv("../results/localSensitivityStatesmTAL.csv")
localSensitivity$X <- NULL
dirLS <- matrix(0, nrow = 40, ncol = 63)

typical <- c(4.75865698e-08, 1.61870719e+02, 9.78576481e-04,
              1.39642352e-03, 3.75000000e-04, 2.12985576e-03, 2.87017424e-03,
              4.28306136e-03, 5.80016806e-05, 1.15774482e-03, 1.37071836e-10,
              3.00666241e-03, 8.55430567e-04, 1.66669761e-05, 6.15845189e-06,
              1.49452825e-07, 3.04893148e-06, 9.19609621e-04, 1.09892277e-04,
              3.95021332e-04, 7.98535285e-03, 1.24649041e-07, 1.22041727e-01,
              2.06235574e-03, 6.69000000e-05, 2.14000000e-02, 1.06908406e-03,
              2.45292416e-03, 2.75659464e-04, 2.03708703e-05, 2.43619703e-03,
              6.30960000e-08, 1.00300000e-03, 1.30000000e-01, 2.44615283e-03,
              2.82430794e-04, 2.43795716e-03, 6.30960000e-08, 4.00000000e-04,
              1.10000000e-01, 1.84518413e-04, 1.53145190e-04, 1.53750000e-04,
              1.33418576e-04, 1.33418576e-04, 3.33216527e-06, 3.33216527e-06,
              5.08045526e-04, 5.08045526e-04, 2.13734943e-04, 2.13734943e-04,
              4.03102759e-05, 4.03102759e-05, 5.44158671e-03, 5.44158671e-03,
              1.50000000e-05, 6.25000000e-05, 7.08405956e-05, 7.08405956e-05,
              7.25000000e-03, 9.12500000e-05, 6.27893000e-03, 2.14431407e-05)

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

## Figure 3.3
pdf("../dataVis/localSensitivityOverallmTAL.pdf", width = 20, height = 20)
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
colnames(dirLS) <- c("PDH Activity", "ICITS Activity", "ACON Activity", 
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
barplot(dirLS[35, ]) ## ATP_c
barplot(dirLS[2, ]) ## dPsi
barplot(dirLS[9, ]) ## NADH_x
barplot(dirLS[10, ]) ## QH2_x
barplot(dirLS[27, ]) ## Cred_i

dirLSsmall <- dirLS[c(2, 9, 10, 27, 35), colSums(abs(dirLS[ 
  c(2, 9, 10, 27, 35), ])) > 0.5]

rownames(dirLSsmall) <- c("Electrical Potential", "NADH", "QH2", 
                          "Reduced CytC", "Cytosolic ATP")
#colnames(dirLSsmall) <- c("Complex IV Activity",
#                          "Total CytC", "Total CoQ", "Total NAD+/NADH", 
#                          "k_O2", "K Leak")

pdf("../dataVis/importantLocalSensitivitymTAL.pdf", width = 10)
par(mar = c(8.1, 10.1, 1.1, 4.1), las = 2, cex.axis = 1)
plot(dirLSsmall, 
     col = viridis::viridis(10), main = "",
     xlab = "", ylab = "")
dev.off()


