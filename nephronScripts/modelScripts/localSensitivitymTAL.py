import scipy.integrate as sci
import numpy as np
import pandas as pd
from SALib.sample import saltelli
from SALib.analyze import sobol

import equations
import pc

ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

def f(t, y, pW, J_AtC, glyc, params): ## Differential equations, with optional arguments specified
    #print(t)
    return equations.conservationEqs1(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType,
                              potassiumW = pW,
                              glyc = glyc,
                              params = params)

normpcPC = pc.pcPC
normParam = pc.paramsmTAL
pc.finalConditions[pc.pcIS.iNADH_x] = pc.pcPC.NADtotmTAL/2.0


def modelPredictions(params):
    ran = np.concatenate((np.linspace(0, 14, 15),
           np.array([16, 17]), np.array([20]),
           np.array([24]), np.linspace(30, 40, 11)), dtype = np.int64,
           casting = 'unsafe')
    ## Change pc.params as necessary
    newParams = np.array(pc.paramsmTAL)
    for i in range(len(ran)):
        newParams[ran[i]] = params[i]
    ## Change pcPC
    pc.pcPC.W_m = params[30]
    pc.pcPC.Ctot = params[31]
    pc.pcPC.Qtot = params[32]
    pc.pcPC.NADtotmTAL = params[33]
    pc.pcPC.FADtot = params[34]
    pc.pcPC.k_O2mTAL = params[35]
    pc.pcPC.k_mADP = params[36]
    ## Change conservationEqs1 arguments
    potassiumWeight = params[37]
    J_AtC = params[38]
    glyc = params[39]

    ## Prepare function
    g = lambda t, y: f(t, y, potassiumWeight, J_AtC, glyc, newParams)

    finalConditions = [4.19130030e-08, 1.65734755e+02, 3.72759315e-04,
                        2.00224069e-03, 3.75000000e-04, 7.98312307e-04, 4.20171769e-03,
                        5.06772854e-04, 7.66939011e-05, 1.36082439e-03, 4.42544075e-11,
                        3.00126707e-03, 6.27337280e-04, 1.13580260e-05, 2.10752403e-06,
                        3.10596564e-07, 8.28312971e-06, 6.05959910e-04, 2.38489939e-04,
                        1.08188813e-03, 8.84143136e-03, 1.35346862e-07, 1.15100483e-01,
                        2.56891026e-03, 6.69000000e-05, 2.14000000e-02, 7.53949892e-04,
                        2.58457683e-03, 1.58847143e-04, 6.41973737e-06, 2.26562763e-04,
                        6.30960000e-08, 1.00300000e-03, 1.30000000e-01, 2.58279256e-03,
                        1.60631412e-04, 2.27026563e-04, 6.30960000e-08, 4.00000000e-04,
                        1.10000000e-01, 2.27987066e-04, 1.53548369e-04, 1.53750000e-04,
                        5.64354505e-05, 5.64354505e-05, 8.80293054e-07, 8.80293054e-07,
                        2.62912402e-04, 2.62912402e-04, 4.51894544e-04, 4.51894544e-04,
                        4.70964418e-05, 4.70964418e-05, 5.15515188e-03, 5.15515188e-03,
                        1.50000000e-05, 6.25000000e-05, 7.08405956e-05, 7.08405956e-05,
                        7.25000000e-03, 9.12500000e-05, 1.03527877e-02, 6.56930339e-06]

    ## Run model to completion
    results = sci.solve_ivp(fun = g,
                            t_span = (0, 500),
                            y0 = finalConditions,
                            method = "LSODA",
                            atol = 1e-10,
                            rtol = 1e-10)
    results = results.y.transpose()[-1]
    print(results)
    return results

def main():
    ## Get parameter ranges to consider in sensitivity analysis
    ran = np.concatenate((np.linspace(0, 14, 15),  ## We only include parameter values
                          ## that are used.
                          np.array([16, 17]), np.array([20]),
                          np.array([24]), np.linspace(30, 40, 11)),
                         dtype=np.int64, casting="unsafe")
    defaults = np.concatenate((np.array(pc.paramsmTAL)[ran], [pc.pcPC.W_m, pc.pcPC.Ctot,
                                                          pc.pcPC.Qtot, pc.pcPC.NADtotmTAL,
                                                          pc.pcPC.FADtot, pc.pcPC.k_O2/2.0,
                                                          pc.pcPC.k_mADP],
                               [1., 1.70e-3, 0.]))
    tops = 1.005*defaults
    bottoms = 0.995*defaults

    ## Consider these 40*2 cases
    y = np.zeros([len(defaults)*2, len(pc.finalConditions)])
    for i in range(len(defaults)):
        print(i)
        baseline = defaults[i]
        for j in range(2):
            if j == 0:
                defaults[i] = tops[i]
            elif j == 1:
                defaults[i] = bottoms[i]
            y[2*i+j, ] = modelPredictions(defaults)
        defaults[i] = baseline

    ## Record this
    y = pd.DataFrame(y)
    y.to_csv("../results/localSensitivityStatesmTAL.csv")

main()