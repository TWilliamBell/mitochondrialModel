mitDismTAL <- read.csv("../results/tailsMitDismTALSim.csv")
param <- feather::read_feather("../results/iterProd.feather")
mitDismTAL <- mitDismTAL[order(mitDismTAL$nums), ]
param <- param[mitDismTAL$nums, ]

basic <- lm(mitDismTAL$ATP_c ~ param$`0` + param$`1` + param$`2` + param$`3`)
summary(basic)

interact <- lm(mitDismTAL$ATP_c ~ param$`0` + param$`1` + param$`2` + param$`3` +
              param$`0`*param$`1`*param$`2`*param$`3`)
summary(interact)

basic <- lm(mitDismTAL$dPsi ~ param$`0` + param$`1` + param$`2` + param$`3`)
summary(basic)

interact <- lm(mitDismTAL$dPsi ~ param$`0` + param$`1` + param$`2` + param$`3` +
                 param$`0`*param$`1`*param$`2`*param$`3`)
summary(interact)

