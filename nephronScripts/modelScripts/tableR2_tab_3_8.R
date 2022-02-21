## Table 3.8

if (!grepl("nephronScripts/modelScripts", getwd())) {
  setwd("nephronScripts/modelScripts")
}

sink(nullfile())
source("atpResponseMitDis_figs.R")
sink()

r2atp <- summary(LMATP)$r.squared
r2atpmix <- summary(LMATPMix)$r.squared

cat("ATP PT R2", "No interactions:", r2atp, "Full interactions:", r2atpmix, "\n")

r2psi <- summary(LMPSI)$r.squared
r2psimix <- summary(LMPSIMix)$r.squared

cat("dPsi PT R2", "No interactions:", r2psi, "Full interactions:", r2psimix, "\n")

sink(nullfile())
source("atpReponseMitDismTAL_figs.R")
sink()

r2atp <- summary(LMATP)$r.squared
r2atpmix <- summary(LMATPMix)$r.squared

cat("ATP PT R2", "No interactions:", r2atp, "Full interactions:", r2atpmix, "\n")

r2psi <- summary(LMPSI)$r.squared
r2psimix <- summary(LMPSIMix)$r.squared

cat("dPsi PT R2", "No interactions:", r2psi, "Full interactions:", r2psimix)