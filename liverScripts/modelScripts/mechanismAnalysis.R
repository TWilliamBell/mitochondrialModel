if (isFALSE(grepl("liverScripts/modelScripts", getwd()))) {
  setwd("liverScripts/modelscripts")
}

files <- c(paste0('../results/resultsMech', 1:2, '.feather'), 
           "../results/resultsMech3.csv", "../results/resultsMech4.feather",
           "../results/resultsMech5.feather", "../results/resultsMech6.feather",
           "../results/resultsMech7.feather",
           "../results/resultsMech8.feather",
           paste0('../results/resultsMech', 1:3, 'C3.csv'),
           "../results/resultsMech4C3.feather", "../results/resultsMech5C3.feather",
           "../results/resultsMech6C3.feather", "../results/resultsMech7C3.feather",
           "../results/resultsMech8C3.feather")

mechanism <- list(tail(feather::read_feather(files[1]), 1), 
                  tail(feather::read_feather(files[2]), 1),
                  tail(data.table::fread(files[3]), 1),
                  tail(feather::read_feather(files[4]), 1),
                  tail(feather::read_feather(files[5]), 1),
                  tail(feather::read_feather(files[6]), 1),
                  tail(feather::read_feather(files[7]), 1),
                  tail(feather::read_feather(files[8]), 1),
                  tail(data.table::fread(files[9]), 1),
                  tail(data.table::fread(files[10]), 1),
                  tail(data.table::fread(files[11]), 1),
                  tail(feather::read_feather(files[12]), 1),
                  tail(feather::read_feather(files[13]), 1),
                  tail(feather::read_feather(files[14]), 1),
                  tail(feather::read_feather(files[15]), 1),
                  tail(feather::read_feather(files[16]), 1))
mechanism <- sapply(mechanism, function(x) x$ATP_c)

baseline <- tail(read.csv("../results/resultsATP.csv"), 1)$ATP_c
baseline <- c(rep(baseline, 3), 0.0021, rep(baseline, 7), 0.0021, rep(baseline, 4))
comparison <- mechanism/baseline

hypoxia <- round(comparison[1:8], 3)
# PT Hypoxia: 0.19
# J_AtC 0.73, Glyc 0.96, Both 0.61, Both+ATP pool 0.94, Both+OXPHOS 0.39,
# J_AtC+OXPHOS 0.58, Glyc+OXPHOS 0.93

C3 <- round(comparison[9:16], 3)
# PT C3: 0.76
# J_AtC 0.95, Glyc 0.98, Both 0.88, Both+ATP pool 1.12, Both+OXPHOS 0.71,
# J_AtC+OXPHOS 0.86, Glyc+OXPHOS 0.98