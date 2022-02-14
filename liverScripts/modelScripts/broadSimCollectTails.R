if (!grepl("liverScripts/modelScripts", getwd())) {
  setwd("./liverScripts/modelScripts")
}

files <- dir("../results")

files <- files[grepl("BroadATP", files)]

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

write.csv(tailsOrderly, file = "../results/tailsBroadSim.csv")