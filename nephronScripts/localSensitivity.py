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
    return equations.conservationEqs1(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType,
                              potassiumW = pW,
                              glyc = glyc,
                              params = params)

normpcPC = pc.pcPC
normParam = pc.params


def modelPredictions(params):
    ran = np.concatenate((np.linspace(0, 14, 15),
           np.array([16, 17]), np.array([20]),
           np.array([24]), np.linspace(30, 40, 11)), dtype = np.int64,
           casting = 'unsafe')
    ## Change pc.params as necessary
    newParams = np.array(pc.params)
    for i in range(len(ran)):
        newParams[ran[i]] = params[i]
    ## Change pcPC
    pc.pcPC.W_m = params[30]
    pc.pcPC.Ctot = params[31]
    pc.pcPC.Qtot = params[32]
    pc.pcPC.NADtot = params[33]
    pc.pcPC.FADtot = params[34]
    pc.pcPC.k_O2 = params[35]
    pc.pcPC.k_mADP = params[36]
    ## Change conservationEqs1 arguments
    potassiumWeight = params[37]
    J_AtC = params[38]
    glyc = params[39]

    ## Prepare function
    g = lambda t, y: f(t, y, potassiumWeight, J_AtC, glyc, newParams)

    ## Run model to completion
    results = sci.solve_ivp(fun = g,
                            t_span = (0, 500),
                            y0 = pc.finalConditions,
                            method = "LSODA",
                            atol = 1e-9,
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
    defaults = np.concatenate((np.array(pc.params)[ran], [pc.pcPC.W_m, pc.pcPC.Ctot,
                                                          pc.pcPC.Qtot, pc.pcPC.NADtot,
                                                          pc.pcPC.FADtot, pc.pcPC.k_O2,
                                                          pc.pcPC.k_mADP],
                               [1., 1.256e-3, 0.0]))
    tops = 1.005*defaults
    bottoms = 0.995*defaults

    ## Consider these 39*2 cases
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
    y.to_csv("../results/localSensitivityStates.csv")

main()