import scipy.integrate as sci
#import feather
import numpy as np
import pandas as pd
#import time

import fluxes
import pc


J_AtC = 7.0e-4 ## Used by Edwards et al.
glyc = 2.3e-4
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clampe

pc.finalConditions[pc.pcIS.iNADH_x] = pc.pcPC.NADtot/2.0
hleakNorm = pc.params[38]

def f(x, i = 2): ## Differential equations, with optional arguments specified
    pc.params[38] = hleakNorm*i
    return fluxes.fluxes(x, param = pc.params,
                              ExpType = ExpType)

def main():
    for j in range(1, 10+1):
        g = lambda x: f(x, i=j)
        hleak = pd.read_csv("../results/resultsHLeak"+str(j-1)+".csv")
        hleak = np.delete(hleak.tail(1).to_numpy(), [0, 1])
        flxs = g(x = hleak)
        print(j)
        print(flxs[pc.pcFl.J_C4]/2.0)
        print(flxs[pc.pcFl.J_F1])
        print(flxs[pc.pcFl.J_F1]/flxs[pc.pcFl.J_C4])
        print("\n")

main()