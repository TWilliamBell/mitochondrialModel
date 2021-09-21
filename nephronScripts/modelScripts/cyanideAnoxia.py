import scipy.integrate as sci
import feather
import numpy as np
import pandas as pd
#import time

import equations
import pc

ExpType = 1 ## in vivo = Pyruvate in  cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

pc.params[38] = 2.0*pc.params[38]
pc.ics[pc.pcIS.iCred_i] = min(pc.ics[pc.pcIS.iCred_i], pc.pcPC.Ctot/2.0)
#pc.pcPC.Ctot = pc.pcPC.Ctot/2.0
#pc.pcPC.Qtot = pc.pcPC.Qtot/2.0


## Hleak is 10x, C4 activity is reduced to 0.2 of default
def f(t, y): ## Differential equations, with optional arguments specified
    return equations.conservationEqs1(y, J_AtC = 1.256e-3,
                              ExpType = ExpType,
                              StateType = StateType,
                              w = [1.0, 1.0, 0.2, 1.0] ## Cyanide!!!
                              )

def g(t, y): ## Differential equations, with optional arguments specified
    return equations.conservationEqs1(y, J_AtC = 4.393e-4,
                              ExpType = ExpType,
                              StateType = StateType,
                              w = [1.0, 1.0, 0.2, 0.2], ## Cyanide!!!
                              tubule = "mTAL")


def main(): ## Runs differential equation for time span and outputs results to
    # a csv file and a feather file.
    results = sci.solve_ivp(fun = f,
                            t_span = (0, 100000),
                            y0 = pc.ics,
                            method = "LSODA",
                            atol = 1e-8,
                            rtol = 1e-10)
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
    results.to_csv("../results/resultsCyanide.csv")
    feather.write_dataframe(results, "../results/resultsCyanide.feather")

    ## mTAL results
    pc.pcPC.k_O2 = pc.pcPC.k_O2 / 2.0 ## Something of this kind seems capable
    ## of explaining the difference between the mTAL and PT

    pc.ics[pc.pcIS.iNADH_x] = min(pc.ics[pc.pcIS.iNADH_x],
                                  pc.pcPC.NADtotmTAL/2.0)
    results = sci.solve_ivp(fun = g,
                            t_span = (0, 100000),
                            y0 = pc.ics,
                            method = "LSODA",
                            atol = 1e-8,
                            rtol = 1e-10)
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
    results.to_csv("../results/resultsCyanidemTAL.csv")
    feather.write_dataframe(results, "../results/resultsCyanidemTAL.feather")

main()