## In order to run this file, you first need to run hypoxia.py and hleak.py

import numpy as np
import pandas as pd

import pc as pc
import fluxes as fl

J_AtC = 7e-4
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

o2norm = pc.finalConditions[pc.pcIS.iO2_x]
hleaknorm = pc.params[38]
hcnorm = pc.ics[pc.pcIS.iH_c]
kcnorm = pc.ics[pc.pcIS.iK_c]
mgcnorm = pc.ics[pc.pcIS.iMg_c]
rho_m = 1.53e-6
convert = 1e9 * rho_m

pc.ics = pc.finalConditions

def main():
    dictToPD = []

    ## Standard case

    normal = pd.read_csv("../results/resultsATP.csv").tail(1).to_numpy()[0]
    normal = np.delete(normal, [0, 1])
    flVals = fl.fluxes(normal, param=pc.params, ExpType=ExpType, potassiumW=1.0)
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
    dictToPD.append(values)
    print('\n')

    ## The other cases

    for i in range(6):
        pc.params[38] = hleaknorm
        pW = 1.0
        J_AtC = 7e-4

        ## Open files
        if i == 0:
            mpFile = pd.read_csv("../results/resultsHypoxiaExtreme.csv").tail(1).to_numpy()[0]
        elif i == 1:
            mpFile = pd.read_csv("../results/resultsHypoxia0.csv").tail(1).to_numpy()[0]
        elif i == 2:
            mpFile = pd.read_csv("../results/resultsHypoxia4.csv").tail(1).to_numpy()[0]
        elif i == 3:
            mpFile = pd.read_csv("../results/resultsHleak1.csv").tail(1).to_numpy()[0]
        elif i == 4:
            mpFile = pd.read_csv("../results/resultsHleak4.csv").tail(1).to_numpy()[0]
        elif i == 5:
            mpFile = pd.read_csv("../results/resultsHleak9.csv").tail(1).to_numpy()[0]

        mpFile = np.delete(mpFile, [0, 1])

        #print(mpFile[pc.pcIS.idPsi])
        if i == 0:
            mpFile[pc.pcIS.iO2_x] = o2norm/72.
        if i == 1:
            mpFile[pc.pcIS.iO2_x] = o2norm*0.1
        if i == 2:
            mpFile[pc.pcIS.iO2_x] = o2norm*0.5
        if i == 3:
            pc.params[38] = hleaknorm*2
        if i == 4:
            pc.params[38] = hleaknorm*5
        if i == 5:
            pc.params[38] = hleaknorm*10

        ## Calculate fluxes
        flVals = fl.fluxes(mpFile, param = pc.params, ExpType = ExpType,
                           potassiumW = pW)

        ## Collecting relevant results
        JO2 = (flVals[pc.pcFl.J_C4]/2.)*convert
        JATP = flVals[pc.pcFl.J_F1]*convert
        PO = JATP/(2.*JO2)
        ATPc = mpFile[pc.pcIS.iATP_c]
        ADPc = mpFile[pc.pcIS.iADP_c]
        ARatio = ATPc/ADPc
        pmf = mpFile[pc.pcIS.idPsi]+((pc.pcPC.RT*np.log(10))/pc.pcPC.F)*\
              ((np.log10(mpFile[pc.pcIS.iH_i])-np.log10(mpFile[pc.pcIS.iH_x])))

        ## Table values
        values = {"JO2": JO2, "JATP": JATP, "PO": PO, "dPsi": mpFile[pc.pcIS.idPsi],
                  "PMF": pmf, "ATPc": ATPc*1000.,
                  "ATP/ADP Ratio": ARatio}
        for key in values.keys():
            print(np.round(values[key], 2), end = '')
            if not key == "ATP/ADP Ratio":
                print(" & ", end = '')
        dictToPD.append(values)
        print("\n")

    PD = pd.DataFrame.from_dict(dictToPD, orient = "columns")
    #print(PD)

## Produces Table 4.4
main()