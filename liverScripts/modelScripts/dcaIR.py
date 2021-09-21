import scipy.integrate as sci
import feather
import numpy as np
import pandas as pd

import equations
import pc

J_AtC = 7.0e-4
glyc = 2.3e-4
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

normPYR = pc.finalConditions[pc.pcIS.iPYR_c]
normO2 = pc.finalConditions[pc.pcIS.iO2_x]

def f(t, y): ## Differential equations, with optional arguments specified
    print(t)
    return equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = ExpType, DCA = 1.5,
                              StateType = StateType, glyc = glyc)

def g(t, y): ## Differential equations, with optional arguments specified
    print(t)
    return equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType, glyc = glyc)

def main(): ## Runs differential equation for time span and outputs results to
    ## a csv file and a feather file.
    count = 0
    if count == 0:
        final = pd.read_feather("../results/resultsIschemia.feather")
        final = np.delete(final.tail(1).to_numpy(), [0])

    finalIschemia = final


    ## Pool becomes smaller post-ischemia
    pc.pcPC.NADtot = pc.pcPC.NADtot * 0.5
    finalIschemia[pc.pcIS.iATP_c] = finalIschemia[pc.pcIS.iATP_c] * 0.25
    finalIschemia[pc.pcIS.iADP_c] = finalIschemia[pc.pcIS.iADP_c] * 0.25
    finalIschemia[pc.pcIS.iAMP_c] = finalIschemia[pc.pcIS.iAMP_c] * 0.25
    pc.finalConditions = finalIschemia
    pc.finalConditions[pc.pcIS.iO2_x] = normO2

    ## Pool is smaller, single stage
    pc.finalConditions = finalIschemia
    pc.finalConditions[pc.pcIS.iO2_x] = normO2
    pc.finalConditions[pc.pcIS.iPYR_x] = normPYR

    results = sci.solve_ivp(fun=f,
                            t_span=(0, 1000),
                            y0=pc.finalConditions,
                            method="LSODA",
                            atol=1e-8,
                            rtol=1e-8)
    final = results.y.transpose()[-1]
    results = np.concatenate((np.array([results.t]), results.y)).transpose()
    results = pd.DataFrame(results,
                           columns=["t", "H_x", "dPsi", "ATP_x", "ADP_x",
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
                                    "PCr_c", "AMP_c"])
    results.to_csv("../results/resultsReperfusionPoolDCA.csv")

    ## Pyruvate overloading
    pc.finalConditions[pc.pcIS.iPYR_x] = normPYR*3.

    results = sci.solve_ivp(fun=g,
                            t_span=(0, 1000),
                            y0=pc.finalConditions,
                            method="LSODA",
                            atol=1e-8,
                            rtol=1e-8)
    final = results.y.transpose()[-1]
    results = np.concatenate((np.array([results.t]), results.y)).transpose()
    results = pd.DataFrame(results,
                           columns=["t", "H_x", "dPsi", "ATP_x", "ADP_x",
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
                                    "PCr_c", "AMP_c"])
    results.to_csv("../results/resultsReperfusionPoolPYR.csv")


main()