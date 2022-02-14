import scipy.integrate as sci
import feather
import numpy as np
import pandas as pd

import equations
import pc
import fluxes as fl

J_AtC = 7.0e-4
glyc = 2.3e-4
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

normPYR = pc.finalConditions[pc.pcIS.iPYR_c]
normO2 = pc.finalConditions[pc.pcIS.iO2_x]

def f(t, y, w = [1., 1., 1., 1.]): ## Differential equations, with optional
    ## arguments specified
    print(t)
    return equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = ExpType, w = w,
                              StateType = StateType, glyc = glyc)

def g(t, y):
    print(t)
    return equations.conservationEqs(y, J_AtC = J_AtC,
                               ExpType = ExpType, w = [0.85, 0.30,
                                                       0.50, 1.00],
                               StateType = StateType, glyc = glyc)

def h(t, y):
    print(t)
    return equations.conservationEqs(y, J_AtC = J_AtC,
                               ExpType = ExpType, w = [0.925, 0.65,
                                                       0.750, 1.00],
                               StateType = StateType, glyc = glyc)

def main(): ## Runs differential equation for time span and outputs results to
    ## a csv file and a feather file.
    count = 0
    pc.finalConditions[pc.pcIS.iO2_x] = normO2/10.0

    pc.pcPC.NADtot = pc.pcPC.NADtot * 0.5
    pc.finalConditions[pc.pcIS.iATP_c] = pc.finalConditions[pc.pcIS.iATP_c] * 0.25
    pc.finalConditions[pc.pcIS.iADP_c] = pc.finalConditions[pc.pcIS.iADP_c] * 0.25
    pc.finalConditions[pc.pcIS.iAMP_c] = pc.finalConditions[pc.pcIS.iAMP_c] * 0.25

    count = 1
    results = sci.solve_ivp(fun = f,
                            t_span = (0, 500),
                            y0 = pc.finalConditions,
                            method = "LSODA",
                            atol = 1e-10,
                            rtol = 1e-10)

    final = results.y.transpose()[-1]
    results = np.concatenate((np.array([results.t]), results.y)).transpose()
    results = pd.DataFrame(results,
                           columns = ["t", "H_x", "dPsi", "ATP_x", "ADP_x",
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
    results.to_csv("../results/resultsIschemiaLessHypox.csv")

    finalIschemia = final

    ## Pool is smaller, single stage
    pc.finalConditions = finalIschemia
    pc.finalConditions[pc.pcIS.iO2_x] = normO2
    pc.finalConditions[pc.pcIS.iPYR_x] = normPYR

    results = sci.solve_ivp(fun=f,
                            t_span=(0, 500),
                            y0=pc.finalConditions,
                            method="LSODA",
                            atol=1e-10,
                            rtol=1e-10)
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
    results.to_csv("../results/resultsReperfusionPoolLessHypox.csv")

#start = time.time()
main()
#end = time.time()
#print(end-start)

#print(a)
#finalConditions = np.array(a.tail(1)) ## Remove the first element before use
#print(finalConditions)