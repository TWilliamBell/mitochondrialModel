import scipy.integrate as sci
import feather
import numpy as np
import pandas as pd
#import time

import equations
import pc

J_AtC = 1.256e-3 ## Used by Edwards et al.
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

#pc.params[11] = pc.params[11]*0.1
#pc.params[38] = pc.params[38]*1.15

#pc.finalConditions[pc.pcIS.iGLU_c] = pc.finalConditions[pc.pcIS.iGLU_c]

def f(t, y): ## Differential equations, with optional arguments specified
    return equations.conservationEqs1(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType,
                              w = [0.5, 1., 1., 1.],
                              cana = 0.33) #0.33



def main(): ## Runs differential equation for time span and outputs results to
    ## a csv file and a feather file.

    results = sci.solve_ivp(fun = f,
                            t_span = (0, 100000),
                            y0 = pc.finalConditions,
                            method = "LSODA",
                            atol = 1e-8,
                            rtol = 1e-9)
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
    results.to_csv("../results/resultsSGLT2i.csv")
    feather.write_dataframe(results, "../results/resultsSGLT2i.feather")
    return results

a = main()
print(a)

b = np.delete(a.tail(1).to_numpy()[0], 0)
import fluxes as fl
print(fl.fluxes(b,
                param = pc.params,
                ExpType = 1,
                w = [0.5, 1., 1., 1.],
                cana = 0.33)[pc.pcFl.J_gdh])
print(fl.fluxes(b,
                param = pc.params,
                ExpType = 1,
                w = [0.5, 1., 1., 1.],
                cana = 0.33)[pc.pcFl.J_got])
print(fl.fluxes(b,
                param = pc.params,
                ExpType = 1,
                w = [0.5, 1., 1., 1.],
                cana = 0.33)[pc.pcFl.J_ASP_GLU])
print(b[pc.pcIS.iGLU_x])
print(pc.pcPC.NADtot-b[pc.pcIS.iNADH_x])
print(b[pc.pcIS.iH_x])