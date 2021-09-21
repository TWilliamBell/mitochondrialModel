mise::mise()

cleanComplex <- function(path, nchars = 10000000) {
  dat <- readChar(con = path, nchars = nchars)
  dat <- gsub("\\(", "", dat)
  dat <- gsub("j\\)", "", dat)
  dat <- gsub(",", "i,", dat)
  dat <- gsub("\\n", "i
              ", dat)
  dat <- gsub(" ", "", dat)
  dat <- paste0(dat, "i")
  dat
}

eigs <- cleanComplex(path = "./results/perturbedEqui.csv")

eigs <- read.csv(text = eigs,
                 colClasses = "complex",
                 header = F)

eigs <- as.matrix(eigs)

pdf("./dataVis/allEigsPerturb.pdf")

plot(eigs, xlab = "Real Part", 
     ylab = "Imaginary Part")

dev.off()

pdf("./dataVis/nearEigsScatterPerturb.pdf")

for (i in 1:nrow(eigs)) {
  colVal = colors()[(i %% 657)+1]
  if (i == 1) {
    print(plot(eigs[i == 1:nrow(eigs), 
               abs(Re(eigs[i, ])) < 100 & abs(Im(eigs[i, ])) < 100], 
               col = colVal, xlab = "Real Part", 
               ylab = "Imaginary Part", ylim = c(-100, 100),
               xlim = c(-100, 100)))
  }
  else {
    print(points(eigs[i == 1:nrow(eigs), 
               abs(Re(eigs[i, ])) < 100 & abs(Im(eigs[i, ])) < 100], 
               col = colVal, xlab = "Real Part", 
               ylab = "Imaginary Part"))
  }
}

dev.off()

pdf("./dataVis/radiusPerturb.pdf")

hist(log(Mod(eigs)), xlab = "Log of the Eigenvalue Radii",
     probability = T, main = "")

dev.off()

## Density

smallOnes <- eigs[abs(Re(eigs)) < 100 & abs(Im(eigs)) < 100]
convertToR2 <- matrix(data = c(Re(smallOnes), Im(smallOnes)), ncol = 2)
convertToR2 <- convertToR2[complete.cases(convertToR2), ]

window <- spatstat::owin(c(-100, 100), c(-100, 100))

pppo <- spatstat::ppp(x = convertToR2[ , 1],
                      y = convertToR2[ , 2], window)

den <- density(pppo, kernel = "gaussian", 
               edge = T, diggle = F, adjust = 1)

pdf("./dataVis/kernelPerturb.pdf")

par(mar = c(5, 3, 2, 0.5)+0.1)
image(den, main = "", col = viridisLite::viridis(15))
axis(side = 1, at = c(-100, -50, 0, 50, 100), labels = T)
axis(side = 2, at = c(-100, -50, 0, 50, 100), labels = T, pos = -108)
mtext(side = c(1, 2), text = c("Real Part", "Imaginary Part"), 
      line = c(2, 2), at = 0)
par(mar = c(5, 4, 4, 2)+0.1)

dev.off()
