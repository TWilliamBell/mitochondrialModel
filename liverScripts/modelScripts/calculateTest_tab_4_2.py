## In order to run this file, you should first run measureTest.py

import numpy as np
import pandas as pd

import pc
import fluxes

J_AtC = 0.0
ExpType = 2 ## in vitro
StateType = 1

def main():

    rho_m =  1.53e-6 ## Lambert et al. 2003
    convert = 1e9 * rho_m

    ## Data loading and cleaning, based on measureVars.py simulations
    s2 = pd.read_csv("../results/resultsState2KHTest.csv").tail(1).to_numpy()[0]
    s2 = np.delete(s2, [0, 1])

    s3 = pd.read_csv("../results/resultsState3KHTest.csv")
    s3 = s3.to_numpy()

    ## Calculate fluxes
    Js3 = list()

    for i in range(len(s3)):
        s3t = np.delete(s3[i], [0, 1])
        Js3.append(fluxes.fluxes(s3t, param = pc.params, ExpType = ExpType))

    Js2 = fluxes.fluxes(s2, param = pc.params, ExpType = ExpType)
    JO2s2 = Js2[2]
    JATPs2 = Js2[3]

    s3JO2 = [item[2] for item in Js3]
    JO2s3 = max(s3JO2)
    ind3 = np.argmax(s3JO2)
    s3ATP = [item[3] for item in Js3]
    JATPs3 = s3ATP[ind3]

    ## Extract values of interest
    valsS2 = {"JO2" : (JO2s2/2.)*convert, "JATP" : JATPs2*convert}
    valsS3 = {"JO2" : (JO2s3/2.)*convert, "JATP" : JATPs3*convert}

    ## All values used in fitting
    print("The oxygen consumption at maximum in State 2 is:")
    print(valsS2["JO2"])
    print("\n")
    print("The oxygen consumption at maximum in State 3 is:")
    print(valsS3["JO2"])
    print("\n")

main()
