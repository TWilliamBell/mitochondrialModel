## Histograms of the Ratios for the Diseased Case

results <- as.data.frame(feather::read_feather("./results/ratiosDiseasedATP.feather"))

for (i in 1:4) {
  png(paste0("./dataVis/", colnames(results)[i], "ATP.png"))
  hist(results[ , i], probability = T, ylab = "Density",
       xlab = "Ratio Found", main = "")
  dev.off()
}

a <- tibble::tibble("J_C1" = NA, "J_C3" = NA, 
                "J_C4" = NA,  "F0F1-ATPase" = NA)

results <- as.data.frame(rbind(results[1:151, ], a, results[152:237, ], a, a,
                               results[238:253, ]))
iterProd <- as.data.frame(feather::read_feather("./results/iterProd.feather"))

model1 <- lm(as.matrix(results[ , 3]) ~ iterProd[ , 1]+
                      iterProd[ , 1]*iterProd[ , 2]+
                      iterProd[ , 1]*iterProd[ , 3]+
                      iterProd[ , 1]*iterProd[ , 4])

sumMod <- summary(model1)

singleVar <- matrix(nrow = 12, ncol = 4)
condensedIterProd <- singleVar

j <- 1

for (i in 1:256) {
  if (sum(iterProd[i, ] == 1) == 3) {
    singleVar[j, ] <- unlist(results[i, ])
    condensedIterProd[j, ] <- unlist(iterProd[i, ])
    j <- j+1
  }
}

datSingle <- cbind(condensedIterProd, singleVar)
iter <- paste0(colnames(results), " parameter scaling")
colnames(datSingle) <- c(iter, colnames(results))

latexV <- xtable::xtable(datSingle)

library(ggplot2)
library(scales)

newDat <- cbind(iterProd, results)

png(paste0("./dataVis/", colnames(results)[1], "StackedATP.png"))
ggplot(data = newDat, aes(x = J_C1, fill = as.factor(newDat[ , 1]))) + 
  geom_histogram(aes(y = ..density..)) + 
  scale_fill_discrete(name = "C1 Weight") +
  ylab("Density")
dev.off()

png(paste0("./dataVis/", colnames(results)[2], "StackedATP.png"))
ggplot(data = newDat, aes(x = J_C3, fill = as.factor(newDat[ , 2]))) + 
  geom_histogram(aes(y = ..density..)) + 
  scale_fill_discrete(name = "C3 Weight") +
  ylab("Density")
dev.off()

png(paste0("./dataVis/", colnames(results)[3], "StackedATP.png"))
ggplot(data = newDat, aes(x = J_C4, fill = as.factor(newDat[ , 3]))) + 
  geom_histogram(aes(y = ..density..)) + 
  scale_fill_discrete(name = "C4 Weight") +
  ylab("Density")
dev.off()

png(paste0("./dataVis/", colnames(results)[4], "StackedATP.png"))
ggplot(data = newDat, aes(x = newDat[ , 8], 
                          fill = as.factor(newDat[ , 1]))) + 
  geom_histogram(aes(y = ..density..)) +
  scale_fill_discrete(name = "C1 Weight") +
  xlab("F0F1-ATPase") +
  ylab("Density")
dev.off()
