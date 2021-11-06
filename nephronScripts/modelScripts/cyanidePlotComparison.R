cyanidePT <- feather::read_feather("./results/resultsCyanide.feather")
baselinesPT <- feather::read_feather("./results/resultsFe.feather")

cyanidemTAL <- feather::read_feather("./results/resultsCyanidemTAL.feather")
baselinesTAL <- feather::read_feather("./results/resultsmTAL.feather")

shortcyanidePT <- cyanidePT[cyanidePT$t < 200, ]
shortbaselinesPT <- baselinesPT[baselinesPT$t < 200, ]

shortcyanidemTAL <- cyanidemTAL[cyanidemTAL$t < 200, ]
shortbaselinesmTAL <- baselinesTAL[baselinesTAL$t < 200, ]

## dPsi

plot(shortbaselinesPT$t, shortbaselinesPT$dPsi, cex = 0.1, ylim = c(140, 400),
     xlab = "Time (s)", ylab = "dPsi (mV)")
lines(shortbaselinesPT$t, shortbaselinesPT$dPsi)
points(shortcyanidePT$t, shortcyanidePT$dPsi, cex = 0.1, col = "red")
lines(shortcyanidePT$t, shortcyanidePT$dPsi, col = "red")

plot(shortbaselinesmTAL$t, shortbaselinesmTAL$dPsi, cex = 0.1, ylim = c(140, 400),
     xlab = "Time (s)", ylab = "dPsi (mV)")
lines(shortbaselinesmTAL$t, shortbaselinesmTAL$dPsi)
points(shortcyanidemTAL$t, shortcyanidemTAL$dPsi, cex = 0.1, col = "red")
lines(shortcyanidemTAL$t, shortcyanidemTAL$dPsi, col = "red")


## Zoomed in
plot(shortbaselinesPT$t, shortbaselinesPT$dPsi, cex = 0.1, ylim = c(140, 200),
     xlab = "Time (s)", ylab = "dPsi (mV)")
lines(shortbaselinesPT$t, shortbaselinesPT$dPsi)
points(shortcyanidePT$t, shortcyanidePT$dPsi, cex = 0.1, col = "red")
lines(shortcyanidePT$t, shortcyanidePT$dPsi, col = "red")

plot(shortbaselinesmTAL$t, shortbaselinesmTAL$dPsi, cex = 0.1, ylim = c(140, 200),
     xlab = "Time (s)", ylab = "dPsi (mV)")
lines(shortbaselinesmTAL$t, shortbaselinesmTAL$dPsi)
points(shortcyanidemTAL$t, shortcyanidemTAL$dPsi, cex = 0.1, col = "red")
lines(shortcyanidemTAL$t, shortcyanidemTAL$dPsi, col = "red")

## Relative

plot(shortbaselinesPT$t, shortbaselinesPT$dPsi/160.2067, cex = 0.1, ylim = c(0.5, 1.5),
     xlab = "Time (s)", ylab = "dPsi (relative to baseline)")
lines(shortbaselinesPT$t, shortbaselinesPT$dPsi/160.2067)
points(shortcyanidePT$t, shortcyanidePT$dPsi/160.2067, cex = 0.1, col = "red")
lines(shortcyanidePT$t, shortcyanidePT$dPsi/160.2067, col = "red")

plot(shortbaselinesmTAL$t, shortbaselinesmTAL$dPsi/165.541774, cex = 0.1, ylim = c(0.5, 1.5),
     xlab = "Time (s)", ylab = "dPsi (relative to baseline)")
lines(shortbaselinesmTAL$t, shortbaselinesmTAL$dPsi/165.541774)
points(shortcyanidemTAL$t, shortcyanidemTAL$dPsi/165.541774, cex = 0.1, col = "red")
lines(shortcyanidemTAL$t, shortcyanidemTAL$dPsi/165.541774, col = "red")

## ATP

plot(shortbaselinesPT$t, shortbaselinesPT$ATP_c, cex = 0.1, 
     ylim = c(0.000, 0.0026),
     xlab = "Time (s)", ylab = "ATP_c (mM)")
lines(shortbaselinesPT$t, shortbaselinesPT$ATP_c)
points(shortcyanidePT$t, shortcyanidePT$ATP_c, cex = 0.1, col = "red")
lines(shortcyanidePT$t, shortcyanidePT$ATP_c, col = "red")

plot(shortbaselinesmTAL$t, shortbaselinesmTAL$ATP_c, cex = 0.1, 
     ylim = c(0.000, 0.0026),
     xlab = "Time (s)", ylab = "ATP_c (mM)")
lines(shortbaselinesmTAL$t, shortbaselinesmTAL$ATP_c)
points(shortcyanidemTAL$t, shortcyanidemTAL$ATP_c, cex = 0.1, col = "red")
lines(shortcyanidemTAL$t, shortcyanidemTAL$ATP_c, col = "red")

## NADH

plot(shortbaselinesPT$t, shortbaselinesPT$NADH_x/0.824e-3, cex = 0.1, ylim = c(0, 1),
     xlab = "Time (s)", ylab = "Redox Ratio")
lines(shortbaselinesPT$t, shortbaselinesPT$NADH_x/0.824e-3)
points(shortcyanidePT$t, shortcyanidePT$NADH_x/0.824e-3, cex = 0.1, col = "red")
lines(shortcyanidePT$t, shortcyanidePT$NADH_x/0.824e-3, col = "red")

plot(shortbaselinesmTAL$t, shortbaselinesmTAL$NADH_x/0.824e-4, cex = 0.1, 
     ylim = c(0, 1),
     xlab = "Time (s)", ylab = "Redox Ratio")
lines(shortbaselinesmTAL$t, shortbaselinesmTAL$NADH_x/0.824e-4)
points(shortcyanidemTAL$t, shortcyanidemTAL$NADH_x/0.824e-4, cex = 0.1, col = "red")
lines(shortcyanidemTAL$t, shortcyanidemTAL$NADH_x/0.824e-4, col = "red")

## Redox of Cytochrome C
plot(shortbaselinesmTAL$t, shortbaselinesmTAL$Cred_i/(1.956e-3), cex = 0.1,
     ylim = c(0, 1), ylab = "Cytochrome C Redox State",
     xlab = "Time (s)")
lines(shortbaselinesmTAL$t, shortbaselinesmTAL$Cred_i/(1.956e-3))
points(shortcyanidemTAL$t, shortcyanidemTAL$Cred_i/(1.956e-3), cex = 0.1,
     col = "red")
lines(shortcyanidemTAL$t, shortcyanidemTAL$Cred_i/(1.956e-3), col = "red")

## Redox of Q protein
plot(shortbaselinesmTAL$t, shortbaselinesmTAL$QH2_x/(2.148e-3), cex = 0.1,
     ylim = c(0, 1), ylab = "Q Protein Redox State",
     xlab = "Time (s)")
lines(shortbaselinesmTAL$t, shortbaselinesmTAL$QH2_x/(2.148e-3))
points(shortcyanidemTAL$t, shortcyanidemTAL$QH2_x/(2.148e-3), cex = 0.1, 
       col = "red")
lines(shortcyanidemTAL$t, shortcyanidemTAL$QH2_x/(2.148e-3), cex = 0.1, 
      col = "red")

## Unsupervised difference searching

PTNorm <- unlist(tail(baselinesPT, 1))
PTCyan <- unlist(tail(cyanidePT, 1))

mTALNorm <- unlist(tail(baselinesTAL, 1))
mTALCyan <- unlist(tail(cyanidemTAL, 1))

PTCyan/PTNorm
mTALCyan/mTALNorm
