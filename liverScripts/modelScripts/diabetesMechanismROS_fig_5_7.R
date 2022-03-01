## Before running this file, run diabetesComplexMechanism.py

if (isFALSE(grepl("liverScripts/modelscripts", getwd()))) {
  setwd("liverScripts/modelscripts")
}

ks <- 1:450

files <- paste0("../results/resultsDiabetesComplexMechanism", ks, ".csv")
## Also compare resultsDiabetesComplexMechanism files from the script diabetesComplexMechanism2.py
files <- lapply(files, data.table::fread)

a <- function(x) tail(x$ATP_c, 1)
q <- function(x) tail(x$QH2_x, 1)
cr <- function(x) tail(x$Cred_i, 1)

ATPs <- sapply(files, a)
Qs <- sapply(files, q)/(2.148e-3*0.2)
Cs <- sapply(files, cr)/(1.956e-3*1.33)

## Figure 5.7
pdf("../dataVis/mechanismDiabetes.pdf", width = 17)
par(mfrow = c(1, 3), cex.lab = 2, cex.axis = 2, mar = c(5.1, 5, 4.1, 2.1))
hist(ATPs*1000, main = "", xlab = "Cytosolic ATP (mM)", xlim = c(0, 8))
hist(Qs, main = "", xlim = c(0, 1), xlab = "Reduced Proportion of Coenzyme Q")
hist(Cs, main = "", xlim = c(0, 1), xlab = "Reduced Proportion of Cytochrome C")
dev.off()
## Give us a sense of the liver's response to these stimuli if you used CytC/CoQ values taken from the
## kidney, which indicates that these values explain the difference in the kidney's response.
