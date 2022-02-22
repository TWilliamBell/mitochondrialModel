## In order to run this file, should first run the script ischemia.py and ischemiamTAL.py

import numpy as np
import pandas as pd

import pc as pc
import fluxes as fl
import fluxesmTAL as flt

J_AtC = 1.256e-3 ## Used by Edwards et al.
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

o2norm = pc.ics[pc.pcIS.iO2_x]
hleaknorm = pc.params[38]
hcnorm = pc.ics[pc.pcIS.iH_c]
kcnorm = pc.ics[pc.pcIS.iK_c]
mgcnorm = pc.ics[pc.pcIS.iMg_c]

## Conversion factor taken from Edwards et al 2019
rho_m = 3.0545e-6
convert = 1e9 * rho_m

def main():
    dictToPD = []

    for i in range(2):
        pc.ics[pc.pcIS.iO2_x] = o2norm
        pc.ics[pc.pcIS.iH_c] = hcnorm
        pc.ics[pc.pcIS.iK_c] = kcnorm
        pc.ics[pc.pcIS.iMg_c] = mgcnorm
        pc.params[38] = hleaknorm
        pW = 1.0

        ## Open files
        if i == 0:
            mpFile = pd.read_csv("../results/resultsReperfusionPool12.csv")
        if i == 1:
            mpFile = pd.read_csv("../results/resultsReperfusionPool12mTAL.csv")
        mpFile = mpFile.tail(1).to_numpy()[0]
        mpFile = np.delete(mpFile, [0, 1])

        ## Conditions to reproduce (for the altered model) Table 3's observations from Edwards et al
        if i == 0:
            J_AtC = 1.23e-3
        if i == 1:
            J_AtC = 1.7e-3

        ## Calculate fluxes
        if i == 0:
            flVals = fl.fluxes(mpFile, param = pc.params, ExpType = ExpType, potassiumW= pW)
        if i == 1:
            flVals = flt.fluxesmTAL(mpFile, param = pc.params, ExpType = ExpType, potassiumW= pW)

        ## Collecting relevant results
        JO2 = (flVals[pc.pcFl.J_C4]/2.)*convert
        JATP = flVals[pc.pcFl.J_F1]*convert
        PO = JATP/(2.*JO2)
        ATPc = mpFile[pc.pcIS.iATP_c]
        ADPc = mpFile[pc.pcIS.iADP_c]
        ARatio = ATPc/ADPc
        pmf = mpFile[pc.pcIS.idPsi] - ((pc.pcPC.RT * np.log(10)) / pc.pcPC.F) * \
              (np.log10(mpFile[pc.pcIS.iH_x]) - np.log10(mpFile[pc.pcIS.iH_c]))
        if i == 0:
            NAD = pc.pcPC.NADtot
        if i == 1:
            NAD = pc.pcPC.NADtotmTAL
        ## Table values
        values = {"dPsi": mpFile[pc.pcIS.idPsi],
                  "PMF": pmf, "ATP/ADP Ratio": ARatio,
                  "NADHRatio": mpFile[pc.pcIS.iNADH_x]/NAD,
                  "CoQRatio": mpFile[pc.pcIS.iQH2_x]/pc.pcPC.Qtot,
                  "CytCRatio": mpFile[pc.pcIS.iCred_i]/pc.pcPC.Ctot}

        for key in values.keys():
            print(np.round(values[key], 2), end = '')
            print(" & ", end = '')
        print("\n")
        dictToPD.append(values)

    PD = pd.DataFrame.from_dict(dictToPD, orient = "columns")
    print(PD)


main()
