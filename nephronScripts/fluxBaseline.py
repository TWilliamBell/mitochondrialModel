import numpy as np
import pandas as pd

import fluxes as fl
import fluxesmTAL as flTAL
import pc

concentrations = np.array(pd.read_feather("../results/resultsATP.feather")
                          .tail(1))[0]

#print(concentrations)

concentrations = np.delete(concentrations, 0)

PTCoA = concentrations[pc.pcIS.iCOASH_x]
PTACCOA = concentrations[pc.pcIS.iACCOA_x]
PTCIT = concentrations[pc.pcIS.iCIT_x]
PTICIT = concentrations[pc.pcIS.iICIT_x]
PTAKG = concentrations[pc.pcIS.iAKG_x]
PTSCOA = concentrations[pc.pcIS.iSCOA_x]
PTSUC = concentrations[pc.pcIS.iSUC_x]
PTFUM = concentrations[pc.pcIS.iFUM_x]
PTMAL = concentrations[pc.pcIS.iMAL_x]
PTOAA = concentrations[pc.pcIS.iOAA_x]

## Fluxes during mitochondrial disease of ATP Synthase
# concentrationsMitDis = pd.read_csv("../results/tailsMitDisPT.csv")
# concentrationsMitDis = concentrationsMitDis.sort_values(by = ['nums'])
# w = [0.25, 0.5, 0.75, 1]
# for i in range(256):
#     conc = np.array(concentrationsMitDis.iloc[i])
#     conc = np.delete(conc, [0, 1, 2, 3])
#     flx = fl.fluxes(conc, param = pc.params, ExpType = 1)[pc.pcFl.J_F1]*w[i % 4]
#     if conc[pc.pcIS.iATP_c] < 0.9*0.00233:
#         continue
#     print(flx/0.00328)

print(fl.fluxes(concentrations, param = pc.params, ExpType = 1)
      [pc.pcFl.J_mdh])

concentrations = np.array(pd.read_feather("../results/resultsmTAL.feather")
                          .tail(1))[0]
concentrations = np.delete(concentrations, 0)

print(flTAL.fluxesmTAL(concentrations, param = pc.params,
                       ExpType = 1)[pc.pcFl.J_mdh])
print("\n")

print(concentrations[pc.pcIS.iCOASH_x]/PTCoA)
print(concentrations[pc.pcIS.iACCOA_x]/PTACCOA)
print(concentrations[pc.pcIS.iCIT_x]/PTCIT)
print(concentrations[pc.pcIS.iICIT_x]/PTICIT)
print(concentrations[pc.pcIS.iAKG_x]/PTAKG)
print(concentrations[pc.pcIS.iSCOA_x]/PTSCOA)
print(concentrations[pc.pcIS.iSUC_x]/PTSUC)
print(concentrations[pc.pcIS.iFUM_x]/PTFUM)
print(concentrations[pc.pcIS.iMAL_x]/PTMAL)
print(concentrations[pc.pcIS.iOAA_x]/PTOAA)
