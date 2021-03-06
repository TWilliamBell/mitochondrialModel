import scipy.integrate as sci
import scipy.optimize as sco
import numpy as np
import pandas as pd

import pc
import measureTest as mtpt
import measureTAL as mtal
import calculateTest as ctpt
import calculateTAL as ctal

def costFn(vals):
    KH = vals[0]
    S3 = vals[1]
    RCR = vals[2]
    LR = vals[3]
    PO = vals[4]
    RCRmTAL = vals[5]
    POmTAL = vals[6]
    if PO > POmTAL:
        compare = (PO-POmTAL)/0.01
    else:
        compare = 0.
    return ((S3 - 4.08) / 0.6) ** 2 + ((RCR - 8.5) / 0.9) ** 2 + \
           ((LR - 0.37) / 0.06) ** 2 + \
           ((PO - 1.8) / 0.1) ** 2 + ((KH - 1.05) / 0.025) ** 2 + \
           ((RCRmTAL - 10.0) / 1.6) ** 2 + \
           ((POmTAL - 1.92) / 0.13) ** 2 + compare

def optimFn(optParam):
    print(optParam)

    khpt = optParam[0]
    khmTAL = 1.0
    ko2pt = optParam[1]
    ko2mTAL = 0.5
    hleakpt = optParam[2]
    hleakmTAL = 1.0
    kleakpt = optParam[3]
    kleakmTAL = 1.0
    antpt = optParam[4]
    antmTAL = 1.0

    pW = kleakpt

    global pc
    pc.pcPC.k_O2 = 1.2e-4*ko2pt
    ## The in vivo case
    ## Increase Potassium-Hydrogen Antiporter activity by 100x from 'new normal'
    ## and see if we get that consistent 1.05 increase in dPsi

    pc.params[37] = khpt * 317200.0
    pc.params[38] = 347.4 * hleakpt
    pc.params[34] = 0.00675 * antpt

    #print(pc.params[37])
    #print(pc.params[38])
    #print(pc.params[34])

    ## PT measure
    try:
        mtpt.main(pW)
    except:
        return 25.0
    ## PT calculate
    ptVals = ctpt.main()
    print("Done PT calculations.")
    ## mTAL params

    #pW = kleakmTAL

    #pc.pcPC.k_O2 = 1.2e-4*ko2mTAL

    #pc.params[37] = khmTAL * (4.7580e+06 / 15)
    #pc.params[38] = 347.4 * hleakmTAL
    #pc.params[34] = 0.00675 * antmTAL

    ## mTAL measure
    #mtal.main(pW)
    ## mTAL calculate
    mtalVals = (8.92, 1.97) #ctal.main()
    print("Done mTAL calculations.")
    vals = ptVals + mtalVals
    print(vals)
    costVal = costFn(vals)
    print(costVal)
    return costVal

import sys, os

# Disable
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore
def enablePrint():
    sys.stdout = sys.__stdout__

def main():
    cost = 100000
    paramSet = [1., 1., 1.5, 1., 1.]
    for i in range(0, 7):
        for j in range(0, 7):
            for k in range(0, 9):
                for l in range(0, 9):
                    for m in range(0, 7):

                        blockPrint()
                        a = optimFn([0.8+(1.2/6.)*i, 0.8+(1.2/6.)*j,
                                     1.+(1./10.)*k,
                                     1.+(1./10.)*l, 0.7+(1.2/6.)*m])
                        enablePrint()

                        if a < cost:
                            cost = a
                            paramSet = [0.8+(1.2/6.)*i, 0.8+(1.2/6.)*j,
                                     1.+(1./9.)*k,
                                     1.+(1./10.)*l, 0.7+(1.3/9.)*m]
                        print(paramSet)
                        print(a)


main()

#optimFn([1., 1., 1.5, 1., 1.])

#optimFn([1., 1., 1., 0.5, 1., 1., 1., 1., 1., 1.]) ## Current values

#optimFn([1.00035857, 1.00035857, 1.00035857, 0.50053786,
#         1.00035857, 1.00035857, 1.00035857, 1.00035857,
#