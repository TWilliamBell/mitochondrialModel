## Before running this script, run mechanism.py, iterprod.py, and broadSim.py, and then
## broadSimCollectTails.R

if (isFALSE(grepl("liverScripts/modelscripts", getwd()))) {
  setwd("liverScripts/modelscripts")
}

files <- c(paste0('../results/resultsMech', 1:2, '.csv'), 
           "../results/resultsMech3.csv", "../results/resultsMech4.csv",
           "../results/resultsMech5.csv", "../results/resultsMech6.csv",
           "../results/resultsMech7.csv",
           "../results/resultsMech8.csv",
           paste0('../results/resultsMech', 1:3, 'C3.csv'),
           "../results/resultsMech4C3.csv", "../results/resultsMech5C3.csv",
           "../results/resultsMech6C3.csv", "../results/resultsMech7C3.csv",
           "../results/resultsMech8C3.csv")

mechanism <- list(tail(data.table::fread(files[1]), 1), 
                  tail(data.table::fread(files[2]), 1),
                  tail(data.table::fread(files[3]), 1),
                  tail(data.table::fread(files[4]), 1),
                  tail(data.table::fread(files[5]), 1),
                  tail(data.table::fread(files[6]), 1),
                  tail(data.table::fread(files[7]), 1),
                  tail(data.table::fread(files[8]), 1),
                  tail(data.table::fread(files[9]), 1),
                  tail(data.table::fread(files[10]), 1),
                  tail(data.table::fread(files[11]), 1),
                  tail(data.table::fread(files[12]), 1),
                  tail(data.table::fread(files[13]), 1),
                  tail(data.table::fread(files[14]), 1),
                  tail(data.table::fread(files[15]), 1),
                  tail(data.table::fread(files[16]), 1)
                  )
mechanism <- sapply(mechanism, function(x) x$ATP_c)

baseline <- tail(read.csv("../results/resultsATP.csv"), 1)$ATP_c
baseline <- c(rep(baseline, 3), 0.0021, rep(baseline, 7), 0.0021, rep(baseline, 4))
comparison <- mechanism/baseline

hypoxia <- round(comparison[1:8], 3)
# PT Hypoxia: 0.19
# J_AtC 0.73, Glyc 0.96, Both 0.61, Both+ATP pool 0.94, Both+OXPHOS 0.39,
# J_AtC+OXPHOS 0.58, Glyc+OXPHOS 0.93
cat("Table 4.6\n")
cat("Hypoxic typical hepatocyte:\t\t\t\t\t", round(6.76/6.91, 3), "\n")
cat("PT-like ATP consumption:\t\t\t\t\t", hypoxia[1], '\n')
cat("No glycolysis:\t\t\t\t\t\t\t", hypoxia[2], "\n")
cat("Kidney-like OXPHOS activities:\t\t\t\t\t", hypoxia[6], "\n")
cat("PT-like ATP consumption and no glycolysis:\t\t\t", hypoxia[3], "\n")
cat("PT-like ATP consumption and OXPHOS activities:\t\t\t", hypoxia[7], "\n")
cat("With all three:\t\t\t\t\t\t\t", hypoxia[5], "\n")
cat("PT cell:\t\t\t\t\t\t\t", round(1.3/2.3, 3))

C3 <- round(comparison[9:16], 3)
# PT C3: 0.76
# J_AtC 0.90, Glyc 0.97, Both 0.78, Both+ATP pool 1.14, Both+OXPHOS 0.74,
# J_AtC+OXPHOS 0.88, Glyc+OXPHOS 0.96

typicalCase <- data.table::fread("../results/tailsBroadSim.csv")
numsKey <- data.table::fread("../results/iterprod.csv")
typicalCase$V1 <- typicalCase$t <- NULL
numsKey$V1 <- NULL
typicalCase <- dplyr::inner_join(numsKey, typicalCase, by = 'nums')
typicalCase <- typicalCase[typicalCase$Glycolysis == 0.00023 &
                             typicalCase$PO2 == 1.0 &
                             typicalCase$CI == 1.0 &
                             typicalCase$CIII == 0.25 &
                             typicalCase$CIV == 1.0 &
                             typicalCase$ATP == 1.0, ]$ATP_c

cat("\n\n\nTable 4.4\n")
cat("Liver typical Complex III reaction:\t\t\t\t", round(typicalCase/baseline[1], 3), "\n")
cat("PT-like ATP consumption:\t\t\t\t\t", C3[1], '\n')
cat("Without glycolytic activity:\t\t\t\t\t", C3[2], "\n")
cat("With kidney-like OXPHOS activity:\t\t\t\t", C3[6], "\n")
cat("With PT-like ATP consumption and no glycolytic activity:\t", C3[3], "\n")
cat("With PT-like ATP consumption and OXPHOS activity:\t\t", C3[7], "\n")
cat("With PT-like OXPHOS activity and no glycolysis:\t\t\t", C3[8], "\n")
cat("With PT-like ATP consumption, no glycolytic activity, \n and kidney-like OXPHOS activity:\t\t\t\t", 
    C3[5], "\n")
cat("PT Case:\t\t\t\t\t\t\t", round(1.7/2.3, 3)) ## Can confirm from Figure 
