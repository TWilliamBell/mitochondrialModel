if (!grepl("liverScripts/modelScripts", getwd())) {
  setwd("./liverScripts/modelScripts")
}

baseline <- data.table::fread("../results/resultsATP.csv")

cyt1 <- data.table::fread("../results/resultsCytC0.csv")
cyt2 <- data.table::fread("../results/resultsCytC1.csv")

tails <- rbind(tail(baseline, 1), tail(cyt1, 1), tail(cyt2, 1))

#barplot(tails$ATP_c)

## Reduction state of the coenyme Q pool
barplot(tails$QH2_x/6.49e-3, names.arg = c("Typical", "0.5x Cytochrome C", 
                                           "0.25x Cytochrome C"),
        col = "mistyrose3")

## Reduction state of the Cytochrome C pool
origpar <- par()
origpar$cin <- origpar$cra <- origpar$csi <- origpar$cxy <- 
  origpar$din <- origpar$page <- NULL
pdf("../dataVis/permCytC.pdf", width = 12)
par(cex.axis = 2)
barplot((tails$Cred_i/4.39e-3)*c(1, 2, 4), 
        names.arg = c("Typical", "0.5x Cytochrome C", 
                      "0.25x Cytochrome C"),
        col = "mistyrose3")
par(origpar)
dev.off()
#barplot(tails$NADH_x)
