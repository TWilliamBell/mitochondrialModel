## In order to run this script, should first run modelPredictionsmTAL.py and mainTAL.py

import numpy as np
import pandas as pd

import pc as pc
import fluxesmTAL as fl

J_AtC = 1.70e-3
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

o2norm = pc.ics[pc.pcIS.iO2_x]
hleaknorm = pc.paramsmTAL[38]
hcnorm = pc.ics[pc.pcIS.iH_c]
kcnorm = pc.ics[pc.pcIS.iK_c]
mgcnorm = pc.ics[pc.pcIS.iMg_c]

## Conversion factor taken from Edwards et al 2019
rho_m = 3.0545e-6
convert = 1e9 * rho_m

def main():
    dictToPD = []

    ## Standard case

    normal = pd.read_csv("../results/resultsmTAL.csv").tail(1).to_numpy()[0]
    normal = np.delete(normal, [0, 1])

    flVals = fl.fluxesmTAL(normal,
                       param=pc.paramsmTAL,
                       ExpType=ExpType,
                       potassiumW=1.0)
    JO2 = (flVals[pc.pcFl.J_C4] / 2.) * convert
    JATP = flVals[pc.pcFl.J_F1] * convert
    PO = JATP / (2. * JO2)
    ATPc = normal[pc.pcIS.iATP_c]
    ADPc = normal[pc.pcIS.iADP_c]
    ARatio = ATPc / ADPc
    pmf = normal[pc.pcIS.idPsi] + ((pc.pcPC.RT * np.log(10)) / pc.pcPC.F) * \
          ((np.log10(normal[pc.pcIS.iH_i])) - (np.log10(normal[pc.pcIS.iH_x])))
    values = {"JO2": JO2, "JATP": JATP, "PO": PO, "dPsi": normal[pc.pcIS.idPsi],
              "PMF": pmf, "ATPc": ATPc * 1000.,
              "ATP/ADP Ratio": ARatio}
    for key in values.keys():
        print(np.round(values[key], 2), end='')
        print(" & ", end='')
    print("\n")

    dictToPD.append(values)

    ## The other cases

    for i in range(21):
        pc.ics[pc.pcIS.iO2_x] = o2norm
        pc.ics[pc.pcIS.iH_c] = hcnorm
        pc.ics[pc.pcIS.iK_c] = kcnorm
        pc.ics[pc.pcIS.iMg_c] = mgcnorm
        pc.params[38] = hleaknorm
        pW = 1.0
        J_AtC = 4.39e-4

        ## Open files
        mpFile = pd.read_csv("../results/resultsMPmTAL"+str(i)+".csv").tail(1).to_numpy()[0]
        mpFile = np.delete(mpFile, [0, 1])

        ## Conditions to reproduce (for the altered model) Table 3's observations from Edwards et al
        if i == 0:
            pc.ics[pc.pcIS.iO2_x] = o2norm*(20.0/10.0)
        if i == 1:
            pc.ics[pc.pcIS.iO2_x] = o2norm*(7.5/10.0)
        if i == 2:
            pc.ics[pc.pcIS.iO2_x] = o2norm*(5.0/10.0)
        if i == 3:
            pc.ics[pc.pcIS.iO2_x] = o2norm*(2.5/10.0)
        if i == 4:
            J_AtC = 0.5*J_AtC
        if i == 5:
            J_AtC = J_AtC*1.25
        if i == 6:
            J_AtC = J_AtC*1.5
        if i == 7:
            pc.ics[pc.pcIS.iH_c] = 10**-7.4
        if i == 8:
            pc.ics[pc.pcIS.iH_c] = 10**-7.0
        if i == 9:
            pc.ics[pc.pcIS.iH_c] = 10**-6.8
        if i == 10:
            pc.ics[pc.pcIS.iK_c] = 60.0e-3
        if i == 11:
            pc.ics[pc.pcIS.iK_c] = 140.0e-3
        if i == 12:
            pc.ics[pc.pcIS.iMg_c] = 0.2e-3
        if i == 13:
            pc.ics[pc.pcIS.iMg_c] = 0.8e-3
        if i == 14:
            pc.paramsmTAL[38] = 0.0
            pW = 0.0
        if i == 15:
            pc.paramsmTAL[38] = 0.0
        if i == 16:
            pc.paramsmTAL[38] = 0.0
            pW = 10.0
        if i == 17:
            pc.paramsmTAL[38] = hleaknorm*10.0
            pW = 0.0
        if i == 18:
            pc.paramsmTAL[38] = hleaknorm*10.0
        if i == 19:
            pW = 10.0
        if i == 20:
            pc.paramsmTAL[38] = hleaknorm*10.0
            pW = 10.0

        ## Calculate fluxes
        flVals = fl.fluxesmTAL(mpFile, param = pc.paramsmTAL,
                               ExpType = ExpType, potassiumW= pW)

        ## Collecting relevant results
        JO2 = (flVals[pc.pcFl.J_C4]/2.)*convert
        JATP = flVals[pc.pcFl.J_F1]*convert
        PO = JATP/(2.*JO2)
        ATPc = mpFile[pc.pcIS.iATP_c]
        ADPc = mpFile[pc.pcIS.iADP_c]
        ARatio = ATPc/ADPc
        pmf = mpFile[pc.pcIS.idPsi]+((pc.pcPC.RT*np.log(10))/pc.pcPC.F)*\
              ((np.log10(mpFile[pc.pcIS.iH_i]))-(np.log10(mpFile[pc.pcIS.iH_x])))

        ## Table values
        values = {"JO2": JO2, "JATP": JATP, "PO": PO, "dPsi": mpFile[pc.pcIS.idPsi],
                  "PMF": pmf, "ATPc": ATPc*1000.,
                  "ATP/ADP Ratio": ARatio}
        for key in values.keys():
            print(np.round(values[key], 2), end = '')
            print(" & ", end = '')
        print("\n")
        dictToPD.append(values)

    PD = pd.DataFrame.from_dict(dictToPD, orient = "columns")
    #print(PD)


main()
