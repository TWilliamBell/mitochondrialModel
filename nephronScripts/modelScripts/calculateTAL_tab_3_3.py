## In order to run this file, first run measureTAL.py

import numpy as np
import pandas as pd

import pc
import fluxesmTAL as fluxes

J_AtC = 0
ExpType = 2 ## in vitro
StateType = 1

def main():

    ## Conversion factor taken from Edwards et al 2019
    rho_m = 3.0545e-6
    convert = 1e9 * rho_m

    ## Data loading and cleaning, based on measureVars.py simulations
    s2 = pd.read_csv("../results/resultsState2TAL.csv").tail(1).to_numpy()[0]
    s2 = np.delete(s2, [0, 1])

    s3 = pd.read_csv("../results/resultsState3TAL.csv")
    s3 = s3.to_numpy()

    po = pd.read_csv("../results/resultsPOTAL.csv")
    po = po.to_numpy()

    ## Calculate fluxes
    Js3 = list()
    Jpo = list()

    for i in range(len(s3)):
        s3t = np.delete(s3[i], [0, 1])
        Js3.append(fluxes.fluxesmTAL(s3t, param = pc.paramsmTAL, ExpType = ExpType))

    for i in range(len(po)):
       pot = np.delete(po[i], [0, 1])
       Jpo.append(fluxes.fluxesmTAL(pot, param = pc.paramsmTAL, ExpType = ExpType))

    Js2 = fluxes.fluxesmTAL(s2, param = pc.paramsmTAL, ExpType = ExpType)
    JO2s2 = Js2[2]
    JATPs2 = Js2[3]

    s3JO2 = [item[2] for item in Js3]
    JO2s3 = max(s3JO2)
    ind3 = np.argmax(s3JO2)
    s3ATP = [item[3] for item in Js3]
    JATPs3 = s3ATP[ind3]

    poJO2 = [item[2] for item in Jpo]
    JO2po = max(poJO2)
    indpo = np.argmax(poJO2)
    poATP = [item[3] for item in Jpo]
    JATPpo = poATP[indpo]

    ## Extract values of interest
    valsS2 = {"JO2" : (JO2s2/2.)*convert, "JATP" : JATPs2*convert}
    valsS3 = {"JO2" : (JO2s3/2.)*convert, "JATP" : JATPs3*convert}
    valsPO = {"JO2" : (JO2po/2.)*convert, "JATP" : JATPpo*convert}

    ## Find RC and PO ratios
    PO = valsPO["JATP"]/(2.*valsPO["JO2"])
    RCR = valsS3["JO2"]/valsS2["JO2"]

    ## All values used in fitting
    # print("State 2 resp:")
    # print(valsS2["JO2"])
    # print("\n")
    # print("The oxygen consumption at maximum in State 3 is:")
    # print(valsS3["JO2"])
    # print("\n")
    print("The Respiratory Control Ratio is:")
    print(RCR)
    print("\n")
    print("The P/O ratio is:")
    print(PO)
    return (RCR, PO)

if __name__ == "__main__":
    main()
