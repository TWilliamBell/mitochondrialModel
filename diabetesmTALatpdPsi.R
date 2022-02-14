diabetes <- dir("../results")
diabetes <- diabetes[grepl("DiabetesMiceUncoupling", diabetes) & 
                       !grepl("DiabetesMiceUncouplingmTAL", diabetes)]


diabetes <- paste0("../results/", diabetes)
diabetes <- lapply(diabetes, data.table::fread)
diabetes <- lapply(diabetes, function(x) tail(x, 1))
ATP <- sapply(diabetes, function(x) x$ATP_c)

dPsi <- sapply(diabetes, function(x) x$dPsi)

pdf("../dataVis/miceCasesPT.pdf")
hist(ATP*1000, main = "",
     xlab = "Cytosolic ATP (mM)", xlim = c(0, 3))
dev.off()

cols <- c("blue", "red")
a <- rep(c(rep(1, 4), rep(2, 4)), 3)
cols <- cols[a]
b <- c(rep(1, 8), rep(2, 8), rep(3, 8))
plot(dPsi, ATP*1000, ylab = "Cytosolic ATP (mM)", 
     xlab = "Electrical Potential Gradient (mV)", col = cols)

diabetes <- dir("../results")
diabetes <- diabetes[grepl("DiabetesMiceUncouplingmTAL", diabetes)]

diabetes <- paste0("../results/", diabetes)
diabetes <- lapply(diabetes, data.table::fread)
diabetes <- lapply(diabetes, function(x) tail(x, 1))
ATP <- sapply(diabetes, function(x) x$ATP_c)
dPsi <- sapply(diabetes, function(x) x$dPsi)

pdf("../dataVis/miceCasesmTAL.pdf")
hist(ATP*1000, main = "",
     xlab = "Cytosolic ATP (mM)", xlim = c(0, 3))
dev.off()

cols0 <- c("blue", "green", "orange", "red")
cols1 <- c("blue", "green", "red")
cols2 <- c("blue", "red")
a <- rep(c(rep(1, 4), rep(2, 4)), 4)
b <- c(rep(1, 8), rep(2, 8), rep(3, 8))
d <- rep(1:4, 4)
cols0 <- cols0[d]
cols1 <- cols1[b]
cols2 <- cols2[a]
plot(dPsi, ATP*1000, ylab = "Cytosolic ATP (mM)", 
     xlab = "Electrical Potential Gradient (mV)", col = cols2)
