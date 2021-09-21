import scipy.integrate as sci
import numpy as np
import pandas as pd
#import time

import equations
import pc
import fluxes

J_AtC = 1.256e-3 ## Used by Edwards et al.
ExpType = 2 ## in vitro
StateType = 1

def s2(t, y): ## For State 2 Resp.
    print(t)
    a = equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = 2,
                              StateType = 1,
                              w = [1., 1., 1., 0.])
    return a

def s3(t, y): ## Differential equations, with optional arguments specified
    print(t)
    return equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = 2,
                              StateType = 1)

def lr(t, y): ## Differential equations, with optional arguments specified
    print(t)
    a = equations.conservationEqs(y, J_AtC = 0.,
                              ExpType = 2,
                              StateType = 1,
                              w = [1., 1., 1., 0.])
    a[pc.pcIS.iADP_c] = 0
    a[pc.pcIS.iATP_c] = 0
    return a

def po(t, y):
    print(t)
    a = equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = 2,
                              StateType = 1)
    a[pc.pcIS.iADP_c] = 0
    a[pc.pcIS.iATP_c] = 0
    return a


def main(): ## Runs differential equation for time span and outputs results to
    ## a csv file and a feather file.

    originalICs = pc.vitroics

    ## State 2 Respiration

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

    resultsS2 = sci.solve_ivp(fun = s2,
                              t_span = (0, 120),
                              y0 = pc.vitroics,
                              method = "LSODA",
                              atol = 1e-6,
                              rtol = 1e-6)

    # results = sci.solve_ivp(fun = s3,
    #                         t_span = (0, 120),
    #                         y0 = pc.vitroics,
    #                         method = "LSODA",
    #                         atol = 1e-8,
    #                         rtol = 1e-8)
    #
    # resultsLR = sci.solve_ivp(fun = lr,
    #                           t_span = (0, 120),
    #                           y0 = pc.vitroics,
    #                           method = "LSODA",
    #                           atol = 1e-8,
    #                           rtol = 1e-8)
    #
    #
    # resultsPO = sci.solve_ivp(fun = po,
    #                           t_span = (0, 120),
    #                           y0 = pc.vitroics,
    #                           method = "LSODA",
    #                           atol = 1e-8,
    #                           rtol = 1e-8)
    #
    # state = results.y[-1]
    # J1 = fluxes.fluxes(state, param = pc.params, ExpType = ExpType)
    #

    resultsS2 = np.concatenate((np.array([resultsS2.t]), resultsS2.y)).transpose()
    resultsS2 = pd.DataFrame(resultsS2,
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
    resultsS2.to_csv("../results/resultsState2F1.csv")

    # results = np.concatenate((np.array([results.t]), results.y)).transpose()
    # results = pd.DataFrame(results,
    #                        columns = ["t", "H_x", "dPsi", "ATP_x", "ADP_x",
    #                                   "AMP_x", "GTP_x", "GDP_x", "Pi_x", "NADH_x",
    #                                   "QH2_x", "OAA_x", "ACCOA_x", "CIT_x", "ICIT_x",
    #                                   "AKG_x", "SCOA_x", "COASH_x", "SUC_x", "FUM_x",
    #                                   "MAL_x", "GLU_x", "ASP_x", "K_x", "Mg_x",
    #                                   "O2_x", "CO2tot_x", "Cred_i", "ATP_i", "ADP_i",
    #                                   "AMP_i", "Pi_i", "H_i", "Mg_i", "K_i", "ATP_c",
    #                                   "ADP_c", "Pi_c", "H_c", "Mg_c", "K_c", "PYR_x",
    #                                   "PYR_i", "PYR_c", "CIT_i", "CIT_c", "AKG_i",
    #                                   "AKG_c", "SUC_i", "SUC_c", "MAL_i", "MAL_c",
    #                                   "ASP_i", "ASP_c", "GLU_i", "GLU_c", "FUM_i",
    #                                   "FUM_c", "ICIT_i", "ICIT_c", "GLC_c", "G6P_c",
    #                                   "PCr_c", "AMP_c"])
    # results.to_csv("../results/resultsState2.csv")
    #
    # resultsLR = np.concatenate((np.array([resultsLR.t]), resultsLR.y)).transpose()
    # resultsLR = pd.DataFrame(resultsLR,
    #                        columns=["t", "H_x", "dPsi", "ATP_x", "ADP_x",
    #                                 "AMP_x", "GTP_x", "GDP_x", "Pi_x", "NADH_x",
    #                                 "QH2_x", "OAA_x", "ACCOA_x", "CIT_x", "ICIT_x",
    #                                 "AKG_x", "SCOA_x", "COASH_x", "SUC_x", "FUM_x",
    #                                 "MAL_x", "GLU_x", "ASP_x", "K_x", "Mg_x",
    #                                 "O2_x", "CO2tot_x", "Cred_i", "ATP_i", "ADP_i",
    #                                 "AMP_i", "Pi_i", "H_i", "Mg_i", "K_i", "ATP_c",
    #                                 "ADP_c", "Pi_c", "H_c", "Mg_c", "K_c", "PYR_x",
    #                                 "PYR_i", "PYR_c", "CIT_i", "CIT_c", "AKG_i",
    #                                 "AKG_c", "SUC_i", "SUC_c", "MAL_i", "MAL_c",
    #                                 "ASP_i", "ASP_c", "GLU_i", "GLU_c", "FUM_i",
    #                                 "FUM_c", "ICIT_i", "ICIT_c", "GLC_c", "G6P_c",
    #                                 "PCr_c", "AMP_c"])
    # resultsLR.to_csv("../results/resultsState2LR.csv")
    #
    #
    # resultsPO = np.concatenate((np.array([resultsPO.t]), resultsPO.y)).transpose()
    # resultsPO = pd.DataFrame(resultsPO,
    #                        columns=["t", "H_x", "dPsi", "ATP_x", "ADP_x",
    #                                 "AMP_x", "GTP_x", "GDP_x", "Pi_x", "NADH_x",
    #                                 "QH2_x", "OAA_x", "ACCOA_x", "CIT_x", "ICIT_x",
    #                                 "AKG_x", "SCOA_x", "COASH_x", "SUC_x", "FUM_x",
    #                                 "MAL_x", "GLU_x", "ASP_x", "K_x", "Mg_x",
    #                                 "O2_x", "CO2tot_x", "Cred_i", "ATP_i", "ADP_i",
    #                                 "AMP_i", "Pi_i", "H_i", "Mg_i", "K_i", "ATP_c",
    #                                 "ADP_c", "Pi_c", "H_c", "Mg_c", "K_c", "PYR_x",
    #                                 "PYR_i", "PYR_c", "CIT_i", "CIT_c", "AKG_i",
    #                                 "AKG_c", "SUC_i", "SUC_c", "MAL_i", "MAL_c",
    #                                 "ASP_i", "ASP_c", "GLU_i", "GLU_c", "FUM_i",
    #                                 "FUM_c", "ICIT_i", "ICIT_c", "GLC_c", "G6P_c",
    #                                 "PCr_c", "AMP_c"])
    # resultsPO.to_csv("../results/resultsState2PO.csv")
    print("Done State 2 calculations.")

    ## State 3 Respiration
    # s2newIC = np.delete(results.tail(1).to_numpy()[0], [0])
    #
    # pc.vitroics = s2newIC
    # pc.vitroics[pc.pcIS.iADP_c] = 2.5e-3
    #
    # results = sci.solve_ivp(fun = s3,
    #                         t_span = (0, 120),
    #                         y0 = pc.vitroics,
    #                         method = "LSODA",
    #                         atol = 1e-6,
    #                         rtol = 1e-8)
    #
    #
    # state = results.y[-1]
    # J2 = fluxes.fluxes(state, param = pc.params, ExpType = ExpType)
    #
    # results = np.concatenate((np.array([results.t]), results.y)).transpose()
    # results = pd.DataFrame(results,
    #                        columns = ["t", "H_x", "dPsi", "ATP_x", "ADP_x",
    #                                   "AMP_x", "GTP_x", "GDP_x", "Pi_x", "NADH_x",
    #                                   "QH2_x", "OAA_x", "ACCOA_x", "CIT_x", "ICIT_x",
    #                                   "AKG_x", "SCOA_x", "COASH_x", "SUC_x", "FUM_x",
    #                                   "MAL_x", "GLU_x", "ASP_x", "K_x", "Mg_x",
    #                                   "O2_x", "CO2tot_x", "Cred_i", "ATP_i", "ADP_i",
    #                                   "AMP_i", "Pi_i", "H_i", "Mg_i", "K_i", "ATP_c",
    #                                   "ADP_c", "Pi_c", "H_c", "Mg_c", "K_c", "PYR_x",
    #                                   "PYR_i", "PYR_c", "CIT_i", "CIT_c", "AKG_i",
    #                                   "AKG_c", "SUC_i", "SUC_c", "MAL_i", "MAL_c",
    #                                   "ASP_i", "ASP_c", "GLU_i", "GLU_c", "FUM_i",
    #                                   "FUM_c", "ICIT_i", "ICIT_c", "GLC_c", "G6P_c",
    #                                   "PCr_c", "AMP_c"])
    # results.to_csv("../results/resultsState3.csv")
    # print("Done State 3 calculations.")

    ## Leak Respiration

    # s2newICLR = np.delete(resultsLR.tail(1).to_numpy(), [0])
    #
    # pc.vitroics = s2newICLR
    # pc.vitroics[pc.pcIS.iADP_c] = 2.5e-3
    #
    # results = sci.solve_ivp(fun = lr,
    #                         t_span = (0, 120),
    #                         y0 = pc.vitroics,
    #                         method = "LSODA",
    #                         atol = 1e-8,
    #                         rtol = 1e-10)
    #
    #
    # state = results.y[-1]
    # J3 = fluxes.fluxes(state, param = pc.params, ExpType = ExpType)
    #
    # results = np.concatenate((np.array([results.t]), results.y)).transpose()
    # results = pd.DataFrame(results,
    #                        columns = ["t", "H_x", "dPsi", "ATP_x", "ADP_x",
    #                                   "AMP_x", "GTP_x", "GDP_x", "Pi_x", "NADH_x",
    #                                   "QH2_x", "OAA_x", "ACCOA_x", "CIT_x", "ICIT_x",
    #                                   "AKG_x", "SCOA_x", "COASH_x", "SUC_x", "FUM_x",
    #                                   "MAL_x", "GLU_x", "ASP_x", "K_x", "Mg_x",
    #                                   "O2_x", "CO2tot_x", "Cred_i", "ATP_i", "ADP_i",
    #                                   "AMP_i", "Pi_i", "H_i", "Mg_i", "K_i", "ATP_c",
    #                                   "ADP_c", "Pi_c", "H_c", "Mg_c", "K_c", "PYR_x",
    #                                   "PYR_i", "PYR_c", "CIT_i", "CIT_c", "AKG_i",
    #                                   "AKG_c", "SUC_i", "SUC_c", "MAL_i", "MAL_c",
    #                                   "ASP_i", "ASP_c", "GLU_i", "GLU_c", "FUM_i",
    #                                   "FUM_c", "ICIT_i", "ICIT_c", "GLC_c", "G6P_c",
    #                                   "PCr_c", "AMP_c"])
    # results.to_csv("../results/resultsLeakResp.csv")
    # print("Done Leak Respiration calculations.")

    ## P/O ratio at half-maximum respiration
    # s2newICPO = np.delete(resultsPO.tail(1).to_numpy()[0], [0])
    # pc.vitroics = s2newICPO
    # pc.vitroics[pc.pcIS.iADP_c] = 2.5e-3*0.008
    # pc.vitroics[pc.pcIS.iATP_c] = 2.0e-3
    #
    #
    # results = sci.solve_ivp(fun = po,
    #                         t_span = (0, 120),
    #                         y0 = pc.vitroics,
    #                         method = "LSODA",
    #                         atol = 1e-8,
    #                         rtol = 1e-10)
    #
    #
    # state = results.y[-1]
    # J4 = fluxes.fluxes(state, param = pc.params, ExpType = ExpType)
    #
    # results = np.concatenate((np.array([results.t]), results.y)).transpose()
    # results = pd.DataFrame(results,
    #                        columns = ["t", "H_x", "dPsi", "ATP_x", "ADP_x",
    #                                   "AMP_x", "GTP_x", "GDP_x", "Pi_x", "NADH_x",
    #                                   "QH2_x", "OAA_x", "ACCOA_x", "CIT_x", "ICIT_x",
    #                                   "AKG_x", "SCOA_x", "COASH_x", "SUC_x", "FUM_x",
    #                                   "MAL_x", "GLU_x", "ASP_x", "K_x", "Mg_x",
    #                                   "O2_x", "CO2tot_x", "Cred_i", "ATP_i", "ADP_i",
    #                                   "AMP_i", "Pi_i", "H_i", "Mg_i", "K_i", "ATP_c",
    #                                   "ADP_c", "Pi_c", "H_c", "Mg_c", "K_c", "PYR_x",
    #                                   "PYR_i", "PYR_c", "CIT_i", "CIT_c", "AKG_i",
    #                                   "AKG_c", "SUC_i", "SUC_c", "MAL_i", "MAL_c",
    #                                   "ASP_i", "ASP_c", "GLU_i", "GLU_c", "FUM_i",
    #                                   "FUM_c", "ICIT_i", "ICIT_c", "GLC_c", "G6P_c",
    #                                   "PCr_c", "AMP_c"])
    # results.to_csv("../results/resultsPO.csv")
    # print("Done P/O ratio calculations.")

#start = time.time()
main()
#end = time.time()
#print(end-start)

#finalConditions = np.array(a.tail(1)) ## Remove the first element before use
