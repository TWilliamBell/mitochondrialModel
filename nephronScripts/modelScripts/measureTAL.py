import scipy.integrate as sci
import numpy as np
import pandas as pd

import equations
import pc

J_AtC = 1.70e-3
StateType = 1

#pc.pcPC.k_O2 = pc.pcPC.k_O2/2.0
#pc.params[34] = 1.00000001*pc.params[34]

def s2(t, y, pW): ## For State 2 Resp.
    #print(t)
    a = equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = 2,
                              StateType = StateType,
                              w = [1., 1., 1., 0.],
                              tubule = "mTAL", potassiumW = pW)
    return a

def s3(t, y, pW): ## Differential equations, with optional arguments specified
    #print(t)
    return equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = 2,
                              StateType = StateType,
                              tubule = "mTAL",
                              potassiumW = pW)

def po(t, y, pW):
    #print(t)
    a = equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = 2,
                              StateType = StateType,
                              tubule = "mTAL",
                              potassiumW = pW)
    a[pc.pcIS.iADP_c] = 0
    a[pc.pcIS.iATP_c] = 0
    return a

def main(pW = 1.): ## Runs differential equation for time span and outputs results to
    ## a csv file and a feather file.

    ## State 2 Respiration
    print(pW)
    print(pc.params[34])
    print(pc.params[37])
    print(pc.params[38])
    print(pc.pcPC.k_O2mTAL)
    ## These inital conditions are used for measurements in Edwards et al
    ## 2020 based on bath conditions in Schiffer et al 2018
    pc.vitroics[pc.pcIS.iH_c] = 10**-7.1
    pc.vitroics[pc.pcIS.iH_i] = pc.vitroics[pc.pcIS.iH_c]
    pc.vitroics[pc.pcIS.iMg_c] = 3.0e-3
    pc.vitroics[pc.pcIS.iMg_i] = pc.vitroics[pc.pcIS.iMg_c]
    pc.vitroics[pc.pcIS.iK_c] = 90e-3
    pc.vitroics[pc.pcIS.iK_i] = pc.vitroics[pc.pcIS.iK_c]
    pc.vitroics[pc.pcIS.iPi_c] = 10e-3
    pc.vitroics[pc.pcIS.iPYR_c] = 5e-3
    pc.vitroics[pc.pcIS.iMAL_c] = 2e-3

    results = sci.solve_ivp(fun = lambda t, y: s3(t, y, pW),
                            t_span = (0, 20),
                            y0 = pc.vitroics,
                            method = "LSODA",
                            atol = 1e-10,
                            rtol = 1e-10)

    resultsPO = sci.solve_ivp(fun = lambda t, y: po(t, y, pW),
                              t_span = (0, 20),
                              y0 = pc.vitroics,
                              method = "LSODA",
                              atol = 1e-10,
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
    results.to_csv("../results/resultsState2TAL.csv")

    resultsPO = np.concatenate((np.array([resultsPO.t]), resultsPO.y)).transpose()
    resultsPO = pd.DataFrame(resultsPO,
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
    #print("Done State 2 calculations.")

    ## State 3 Respiration
    s2newIC = np.delete(results.tail(1).to_numpy()[0], [0])

    pc.vitroics = s2newIC
    pc.vitroics[pc.pcIS.iADP_c] = 2.5e-3

    results = sci.solve_ivp(fun = lambda t, y: s3(t, y, pW),
                            t_span = (0, 20),
                            y0 = pc.vitroics,
                            method = "LSODA",
                            atol = 1e-10,
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
    results.to_csv("../results/resultsState3TAL.csv")
    #print("Done State 3 calculations.")

    ## P/O ratio at half-maximum respiration
    s2newICPO = np.delete(resultsPO.tail(1).to_numpy()[0], [0])
    pc.vitroics = s2newICPO
    pc.vitroics[pc.pcIS.iADP_c] = 2.5e-3*0.008
    pc.vitroics[pc.pcIS.iATP_c] = 2.0e-3


    results = sci.solve_ivp(fun = lambda t, y: po(t, y, pW),
                            t_span = (0, 20),
                            y0 = pc.vitroics,
                            method = "LSODA",
                            atol = 1e-10,
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
    results.to_csv("../results/resultsPOTAL.csv")
    #print("Done P/O ratio calculations.")

#main()
