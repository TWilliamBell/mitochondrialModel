if (!grepl("mitochondrialModel/modelScripts", getwd())) {
  setwd("./modelScripts")
}

diabetesmTAL <- list()

files <- dir("../results")

files <- files[grepl("resultsDiabetesComplex", files) & !grepl("mTAL", files)]

nums <- gsub("[^\\d]+", "", files, perl=TRUE)
nums[nums == ""] <- "100000"
nums <- as.numeric(nums)
nums[nums == 1e5] <- NA_real_

tails <- list()
files <- paste0("../results/", files)
for (i in files) {
  tails[[i]] <- as.data.frame(tail(data.table::fread(i), 1))
}

tailsOrderly <- as.data.frame(cbind(nums, do.call(rbind, tails)))
tailsOrderly <- tailsOrderly[order(tailsOrderly$nums), ]

write.csv(tailsOrderly, file = "../results/tailsDiabetesComplex.csv")
