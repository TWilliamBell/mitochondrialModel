## Using feather and dplyr

results <- feather::read_feather("./results/resultsATP.feather")

# par(mfrow = c(3, 2), mar = c(0.5, 0.5, 0.5, 0.5))
# 
# for (i in seq_len(ncol(results)-1)) {
#   print(plot(results$t, dplyr::pull(results, i+1), 
#        cex = 0.1))
#   Sys.sleep(0.1)
# }
# 
# rm(i)
# 
# par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

results <- as.data.frame(results)
final <- unlist(tail(results, 1))[2:64]

