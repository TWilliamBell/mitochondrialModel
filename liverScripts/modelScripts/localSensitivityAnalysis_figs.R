## In order to run this script, first run localSensitivity.py

if (!grepl("liverScripts/modelScripts", getwd())) {
  setwd("liverScripts/modelScripts")
}

localSensitivity <- read.csv("../results/localSensitivityStates.csv")
localSensitivity$X <- NULL
dirLS <- matrix(0, nrow = 40, ncol = 63)

typical <- c( 5.66635757e-08,  1.65429541e+02,  3.75420860e-03,  1.27157914e-02,
              3.75000000e-04,  1.16962298e-03,  3.83040702e-03,  1.42447619e-03,
              6.10862490e-05,  7.93295084e-05,  6.67776873e-11,  3.00924392e-03,
              8.13841220e-04,  1.27370874e-05,  4.17201894e-02,  4.06744389e-08,
              5.76204010e-07,  1.55880188e-06,  1.78506098e-05,  1.01470406e-06,
              1.39789900e-10,  0.00000001,     1.12445693e-01, 3.75695934e-03,
              3.18000000e-05,  2.14000000e-02,  6.76796556e-04,  6.81481652e-03,
              4.65293580e-04,  2.19427318e-05,  1.04284027e-03,  6.30960000e-08,
              1.00300000e-03,  1.30000000e-01,  6.81189472e-03,  4.68215376e-04,
              1.04359976e-03,  6.30960000e-08,  4.00000000e-04,  8.80000000e-02,
              1.67335930e-04,  1.53444863e-04,  1.53750000e-04,  7.09920524e-04,
              7.09920525e-04,  4.99842889e-04,  4.99842867e-04,  1.51603878e-05,
              1.51603898e-05,  1.22065167e-08,  1.22065162e-08,  5.12441942e-04,
              5.12441942e-04,  1.20988662e-09,  1.05922976e-10,  1.50000000e-05,
              6.25000000e-05,  7.11112720e-05,  7.11112720e-05,  7.25000000e-03,
              9.12500000e-05,  9.35861547e-03,  2.22287048e-05)

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

## Figure 4.1
pdf("../dataVis/localSensitivityOverallLiver.pdf", width = 20, height = 20)
par(mar = c(10.1, 10.1, 2.1, 4.1), las = 2)
tempDirLS <- dirLS
rownames(tempDirLS) <- c("PDH Activity", "CITS Activity", "ACON Activity", 
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
                     "K Leak", "Max. ATP Consumption", "Glycolytic ATP")
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
plot(t(tempDirLS), col = viridis::viridis(10), xlab = "", ylab = "", main = "")
dev.off()

pdf("../dataVis/localSensitivityOverallQualitative.pdf", width = 20, height = 20)
par(mar = c(10.1, 10.1, 2.1, 4.1), las = 2)
plot(t(sign(dirLS)*((abs(dirLS))^0.1)), col = viridis::viridis(10), xlab = "", ylab = "", main = "")
dev.off()

importantStatedirLS <- dirLS[ , colSums(abs(dirLS)) != 0]
image(abs(importantStatedirLS))
plot(abs(importantStatedirLS), col = viridis::viridis)

importantParamdirLS <- dirLS[rowSums(abs(dirLS)) != 0, ]
image(abs(importantParamdirLS))
plot(importantParamdirLS, col = viridis::viridis)

vimportantStatedirLS <- dirLS[ , colSums(abs(dirLS)) > 1]
image(abs(importantStatedirLS))
plot(vimportantStatedirLS, col = viridis::viridis)

vimportantParamdirLS <- dirLS[rowSums(abs(dirLS)) > 1, ]
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
importantParams <- paramNames[rowSums(abs(dirLS)) > 1]

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
importantStates <- stateNames[colSums(abs(dirLS)) > 1]

par(mar = c(5.1, 4.1, 4.1, 2.1), cex.axis = 0.5)
barplot(dirLS[ , 35]) ## ATP_c
barplot(dirLS[ , 2]) ## dPsi
barplot(dirLS[ , 9]) ## NADH_x
barplot(dirLS[ , 10]) ## QH2_x
barplot(dirLS[ , 27]) ## Cred_i

dirLSsmall <- dirLS[rowSums(abs(dirLS[ , 
                                       c(2, 9, 10, 27, 35)])) > 0.5, 
                    c(2, 9, 10, 27, 35)]
rownames(dirLSsmall) <- c("PDH Activity", "AKGD Activity",
                          "Complex III Activity",
                          "Complex IV Activity", "k_Pi1",
                          "Total CytC", "Total CoQ", "Total NAD/NADH", "k_O2",
                          "Max. ATP Consumption")
colnames(dirLSsmall) <- c("Electrical Potential", "NADH", "QH2", 
                          "Reduced CytC", "Cytosolic ATP")

## Figure 4.2
pdf("../dataVis/importantLocalSensitivityLiver.pdf", width = 10)
par(mar = c(10.1, 10.1, 1.1, 4.1), las = 2, cex.axis = 1)
plot(t(dirLSsmall), 
     col = viridis::viridis(10), main = "",
     xlab = "", ylab = "")
dev.off()

