## Saved you a run:
# S2 JO2: 4.22, RCR: 9.97, LR: 0.577, PO: 2.12

import numpy as np
import pandas as pd

import pc
import fluxes

J_AtC = 1.256e-3 ## Used by Edwards et al.
ExpType = 2 ## in vitro
StateType = 1

def main():

    ## Conversion factor taken from Edwards et al 2019
    rho_m = 3.0545e-6
    convert = 1e9 * rho_m

    ## Data loading and cleaning, based on measureVars.py simulations
    s2 = pd.read_csv("../results/resultsState2.csv").tail(1).to_numpy()[0]
    s2 = np.delete(s2, [0, 1])
    #s2LR = pd.read_csv("../results/resultsState2LR.csv").tail(1).to_numpy()[0]
    #s2LR = np.delete(s2LR, [0, 1]


    print("State 2 State Variables:")
    print(s2)
    print("\n")

    s3 = pd.read_csv("../results/resultsState3.csv")

    print("State 3 State Variables:")
    print(np.delete(s3.tail(1).to_numpy()[0], [0, 1]))
    print("\n")

    s3 = s3.to_numpy()

    lr = pd.read_csv("../results/resultsLeakResp.csv")

    print("Leak respiration State Variables:")
    print(np.delete(lr.tail(1).to_numpy()[0], [0, 1]))
    print("\n")

    lr = lr.to_numpy()

    po = pd.read_csv("../results/resultsPO.csv")

    print("P/O measurement State Variables:")
    print(np.delete(po.tail(1).to_numpy()[0], [0, 1]))
    print("\n")

    po = po.to_numpy()

    ## Calculate fluxes
    #Js2 = list()
    Js3 = list()
    Jlr = list()
    Jpo = list()

    #for i in range(len(s2)):
    #    s2t = np.delete(s2[i], [0, 1])
    #    Js2.append(fluxes.fluxes(s2t, param = pc.params, ExpType = ExpType))

    for i in range(len(s3)):
        s3t = np.delete(s3[i], [0, 1])
        Js3.append(fluxes.fluxes(s3t, param = pc.params, ExpType = ExpType))

    for i in range(len(lr)):
        lrt = np.delete(lr[i], [0, 1])
        Jlr.append(fluxes.fluxes(lrt, param = pc.params, ExpType = ExpType))

    for i in range(len(po)):
        pot = np.delete(po[i], [0, 1])
        Jpo.append(fluxes.fluxes(pot, param = pc.params, ExpType = ExpType))

    #s2JO2 = [item[2] for item in Js2]
    #JO2s2 = max(s2JO2)
    #ind2 = np.argmax(s2JO2)
    #s2ATP = [item[3] for item in Js2]
    #JATPs2 = s2ATP[ind2]
    Js2 = fluxes.fluxes(s2, param = pc.params, ExpType = ExpType)
    JO2s2 = Js2[2]
    JATPs2 = Js2[3]

    #Js2LR = fluxes.fluxes(s2LR, param = pc.params, ExpType = ExpType)
    #JO2s2LR = Js2LR[2]
    #print("Leak respiration state 2 oxygen consumption:")
    #print(JO2s2LR)

    s3JO2 = [item[2] for item in Js3]
    JO2s3 = max(s3JO2)
    ind3 = np.argmax(s3JO2)
    s3ATP = [item[3] for item in Js3]
    JATPs3 = s3ATP[ind3]


    lrJO2 = [item[2] for item in Jlr]
    JO2lr = max(lrJO2)
    indlr = np.argmax(lrJO2)
    lrATP = [item[3] for item in Jlr]
    JATPlr = lrATP[indlr]

    poJO2 = [item[2] for item in Jpo]
    JO2po = max(poJO2)
    indpo = np.argmax(poJO2)
    poATP = [item[3] for item in Jpo]
    JATPpo = poATP[indpo]



    ## Extract values of interest
    valsS2 = {"JO2" : (JO2s2/2.)*convert, "JATP" : JATPs2*convert}
    valsS3 = {"JO2" : (JO2s3/2.)*convert, "JATP" : JATPs3*convert}
    valsLR = {"JO2" : (JO2lr/2.)*convert, "JATP" : JATPlr*convert}
    valsPO = {"JO2" : (JO2po/2.)*convert, "JATP" : JATPpo*convert}

    ## Find RC and PO ratios
    PO = valsPO["JATP"]/(2*valsPO["JO2"])
    RCR = valsS3["JO2"]/valsS2["JO2"]

    ## All values used in fitting
    print("The oxygen consumption at maximum in State 3 is:")
    print(valsS3["JO2"])
    print("\n")
    print("The Respiratory Control Ratio is:")
    print(RCR)
    print("\n")
    print("The leak respiration oxygen consumption is:")
    print(valsLR["JO2"])
    print("\n")
    print("The P/O ratio is:")
    print(PO)

main()
