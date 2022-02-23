## Table 3.8
## In order to run this file, you should first run simulateMitDis1.py and then collectMitDisPT.R &
## writeProduct.py, then run simulateMitDismTAL.py, and then collectMitDismTAL.R

if (!grepl("nephronScripts/modelScripts", getwd())) {
  setwd("nephronScripts/modelScripts")
}

sink(nullfile())
source("atpResponseMitDis_figs.R")
sink()

r2atp <- summary(LMATP)$r.squared
r2atpmix <- summary(LMATPMix)$r.squared

cat(paste0("ATP PT R2\t", "No interactions:\t", r2atp, "\t Full interactions:\t", r2atpmix, "\n"))

r2psi <- summary(LMPSI)$r.squared
r2psimix <- summary(LMPSIMix)$r.squared

cat(paste0("dPsi PT R2\t", "No interactions:\t", r2psi, "\t Full interactions:\t", r2psimix, "\n"))

sink(nullfile())
suppressMessages(source("atpReponseMitDismTAL_figs.R"))
sink()

r2atp <- summary(LMATP)$r.squared
r2atpmix <- summary(LMATPMix)$r.squared

cat(paste0("ATP mTAL R2\t", "No interactions:\t", r2atp, "\t Full interactions:\t", r2atpmix, "\n"))

r2psi <- summary(LMPSI)$r.squared
r2psimix <- summary(LMPSIMix)$r.squared

cat(paste0("dPsi mTAL R2\t", "No interactions:\t", r2psi, "\t Full interactions:\t", r2psimix))

