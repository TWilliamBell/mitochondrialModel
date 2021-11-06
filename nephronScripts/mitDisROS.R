if (!grepl("mitochondrialModel/modelScripts", getwd())) {
  setwd("./modelScripts")
}

mitDis <- list()

trueVec <- rep(TRUE, 256)

for (i in 1:256) {
  if (paste0("resultsDiseaseATP", i, "Fe.feather") %in% dir("../results")) {
    mitDis[[i]] <- feather::read_feather(paste0("../results/resultsDiseaseATP", 
                                                i, "Fe.feather"))
  } else {
    trueVec[i] <- FALSE
  }
}

normalTail <- tail(read.csv("../results/resultsATP.csv"), 1)

mitDis <- do.call(rbind, lapply(mitDis, function(x) as.data.frame(tail(x, 1))))

Q <- c(mitDis$QH2_x, 0.0002288496)/2.418e-3
C <- c(mitDis$Cred_i, 0.000708212)/1.956e-3
N <- c(mitDis$NADH_x, 6.050095e-05)/0.824e-3

hist(N, main = "", xlab = "Reduced Proportion of NADH/NAD+ Pool",
     xlim = c(0, 1), col = "peachpuff1")

hist(Q, main = "", xlab = "Reduced Proportion of Coenzyme Q Pool",
     xlim = c(0, 1), col = "peachpuff1")
abline(v = normalTail[12]/2.418e-3, col = "red")

hist(C, main = "", xlab = "Reduced Proportion of Cytochrome C Pool",
     xlim = c(0, 1), col = "peachpuff1")
abline(v = normalTail[29]/1.956e-3, col = "red")

