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

def f(t, y, pW, J_AtC, glyc, params):
    ## Differential equations, with optional arguments specified
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
    newParams[ran] = params[0:len(ran)]
    #pc.params = list(newParams)
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
    print(results, flush = True)
    return results

stateVars = ["H_x", "dPsi", "ATP_x", "ADP_x",
            "AMP_x", "GTP_x", "GDP_x", "Pi_x", "NADH_x",
            "QH2_x", "OAA_x", "ACCOA_x", "CIT_x", "ICIT_x",
            "AKG_x", "SCOA_x", "COASH_x", "SUC_x", "FUM_x",
            "MAL_x", "GLU_x", "ASP_x", "K_x", "Mg_x",
            "O2_x", "CO2tot_x", "Cred_i", "ATP_i", "ADP_i",
            "AMP_i", "Pi_i", "H_i", "Mg_i", "K_i", "ATP_c",
            "ADP_c", "Pi_c", "H_c", "Mg_c", "K_c", "PYR_x",
            "PYR_i", "PYR_c", "CIT_i", "CIT_c", "AKG_i",
            "AKG_c", "SUC_i", "SUC_c", "MAL_i", "MAL_c",
            "ASP_i", "ASP_c", "GLU_i", "GLU_c", "FUM_i",
            "FUM_c", "ICIT_i", "ICIT_c", "GLC_c", "G6P_c",
            "PCr_c", "AMP_c"]

unimportantStateVars = ['AMP_x', 'CO2tot_x', 'H_i', 'Mg_i', 'K_i', 'H_c',
                        'Mg_c', 'K_c', 'PYR_c', 'FUM_i', 'FUM_c', 'ICIT_i',
                        'ICIT_c', 'GLC_c', 'G6P_c']

def main():
    ## Make parameter ranges
    print("Starting step 1 of 4", flush = True)
    ran = np.concatenate((np.linspace(0, 14, 15), ## We only include parameter values
                          ## that are used.
                          np.array([16, 17]), np.array([20]),
                          np.array([24]), np.linspace(30, 40, 11)),
                         dtype = np.int64, casting = "unsafe")
    defaults = np.concatenate((np.array(pc.params)[ran], [pc.pcPC.W_m, pc.pcPC.Ctot,
                                                          pc.pcPC.Qtot, pc.pcPC.NADtot,
                                                          pc.pcPC.FADtot, pc.pcPC.k_O2,
                                                          pc.pcPC.k_mADP],
                                                        [1., 1.256e-3, 2.3e-4]))
    bottoms = defaults*0.4
    tops = defaults*2.5
    ranges = [[bottoms[i], tops[i]] for i in range(len(defaults))]

    ## Sample ranges
    print("Starting step 2 of 4", flush = True)
    ranges = {
        'num_vars' : len(defaults),
        'names' : ['Vmax1', 'Vmax2', 'Vmax3', 'Vmax4', 'Vmax5', 'Vmax6',
                   'Vmax7', 'Vmax8', 'Vmax9', 'Vmax10', 'Vmax11', 'x_PYR_H',
                   'x_GLU_H', 'x_CIT_MAL', 'x_AKG_MAL', 'x_MAL_PI',
                   'x_ASP_GLU', 'x_SUC_MAL', 'x_KH', 'x_C1', 'x_C3', 'x_C4',
                   'x_F1', 'x_ANT', 'x_Pi1', 'k_PiH', 'x_KH', 'x_Hle',
                   'k_Pi1', 'k_Pi2', 'W_m', 'Ctot', 'Qtot', 'NADtot',
                   'FADtot', 'k_O2', 'k_mADP', 'pW', 'J_AtC', 'glyc'],
        'bounds' : ranges
    }
    param_values = saltelli.sample(ranges, 256)

    ## Make predictions for samples
    print("Starting step 3 of 4", flush = True)
    y = np.zeros([param_values.shape[0], len(pc.finalConditions)])
    parVals = pd.DataFrame(param_values)
    parVals.to_csv("../results/parameterStatesSmall.csv")
    for i, par in enumerate(param_values):
        print(i, flush = True)
        y[i,:] = modelPredictions(par)
    state = pd.DataFrame(y)
    state.to_csv("../results/sensitivityStatesSmall.csv")

    # state = pd.read_csv("../results/sensitivityStates.csv")
    # y = np.array(state)[:, 1:64]

    ## No missing values
    #anyNANY = np.any([np.isnan(np.sum(y[i, ])) for i in range(y.shape[0])])

    ## Analyze samples
    print("Starting step 4 of 4", flush = True)
    importantStateVars = list()
    stds = y.std(axis = 0)
    for i in range(len(pc.finalConditions)):
        if stds[i] < 1e-10:
            continue
        importantStateVars.append(stateVars[i])
        Si = sobol.analyze(ranges, y[:,i])
        totalSi, firstSi, secondSi = Si.to_df()
        totalSi.to_csv("../results/totalSiSmall"+str(i)+".csv")
        firstSi.to_csv("../results/firstSiSmall"+str(i)+".csv")
        secondSi.to_csv("../results/secondSiSmall"+str(i)+".csv")
    print(importantStateVars)
    print("Done.")

main()
