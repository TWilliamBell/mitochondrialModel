import scipy.integrate as sci
import numpy as np
import pandas as pd
import time

import equations
import pc
import fluxes

J_AtC = 0.

ExpType = 1 ## in vivo
StateType = 1

def kh(t, y, pW):
    return equations.conservationEqs1(y, J_AtC = 1.256e-3,
                              ExpType = ExpType,
                              StateType = StateType,
                              w = [1., 1., 1., 1.],
                              potassiumW = pW)

ExpType = 2 ## Change to in vitro

def s2(t, y, pW): ## For State 2 Resp.
    #print(t)
    a = equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType,
                              w = [1., 1., 1., 0.],
                              potassiumW = pW, timeStart = time.time())
    return a

def s3(t, y, pW): ## Differential equations, with optional arguments specified
    #print(t)
    return equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType,
                              potassiumW = pW, timeStart = time.time())

def lr(t, y, pW): ## Differential equations, with optional arguments specified
    #print(t)
    a = equations.conservationEqs(y, J_AtC = 0.,
                              ExpType = ExpType,
                              StateType = StateType,
                              w = [1., 1., 1., 0.],
                              potassiumW = pW, timeStart = time.time())
    a[pc.pcIS.iADP_c] = 0
    a[pc.pcIS.iATP_c] = 0
    return a

def po(t, y, pW):
    #print(t)
    a = equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType,
                              w = [1., 1., 1., 1.],
                              potassiumW = pW, timeStart = time.time())
    a[pc.pcIS.iADP_c] = 0
    a[pc.pcIS.iATP_c] = 0
    return a

def main(pW = 1.): ## Runs differential equation for time span and outputs results to
    ## a csv file and a feather file.

    ## The in vivo case
    ## Increase Potassium-Hydrogen Antiporter activity by 100x from 'new normal'
    ## and see if we get that consistent 1.05 increase in dPsi

    pc.params[37] = pc.params[37] * 20


    resultsKH = sci.solve_ivp(fun = lambda t, y: kh(t, y, pW = pW),
                              t_span = (0, 20),
                              y0 = pc.finalConditions,
                              method = "Radau", t_eval = np.linspace(0, 20, 100000),# "LSODA",
                              atol = 1e-9,
                              rtol = 1e-9)

    resultsKH = np.concatenate((np.array([resultsKH.t]), resultsKH.y)).transpose()
    resultsKH = pd.DataFrame(resultsKH,
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
    resultsKH.to_csv("../results/resultsKH100Test.csv")

    ## For everything after this point use the following:
    pc.params[37] = pc.params[37] / 20.0

    resultsKHbase = sci.solve_ivp(fun = lambda t, y: kh(t, y, pW = pW),
                              t_span = (0, 20),
                              y0 = pc.finalConditions,
                              method = "Radau", t_eval = np.linspace(0, 20, 100000),# "LSODA",
                              atol = 1e-9,
                              rtol = 1e-9)

    resultsKHbase = np.concatenate((np.array([resultsKHbase.t]), resultsKHbase.y)).transpose()
    resultsKHbase = pd.DataFrame(resultsKHbase,
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
    resultsKHbase.to_csv("../results/resultsKHbaseTest.csv")

    originalICs = pc.vitroics

    ## What follows is for the in vitro case

    ## State 2 Respiration

    print(pW)
    print(pc.params[34])
    print(pc.params[37])
    print(pc.params[38])
    print(pc.pcPC.k_O2)

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

    print("Everything before S2")

    # resultsS2 = sci.solve_ivp(fun = lambda t, y: s2(t, y, pW = pW),
    #                           t_span = (0, 20),
    #                           y0 = pc.vitroics,
    #                           method = "LSODA", t_eval = np.linspace(0, 20, 100000),# "LSODA",
    #                           atol = 1e-10,
    #                           rtol = 1e-10)

    results = sci.solve_ivp(fun = lambda t, y: s3(t, y, pW = pW),
                            t_span = (0, 20),
                            y0 = pc.vitroics,
                            method = "LSODA", t_eval = np.linspace(0, 20, 100000),# "LSODA",
                            atol = 1e-10,
                            rtol = 1e-10)

    resultsLR = sci.solve_ivp(fun = lambda t, y: lr(t, y, pW = pW),
                              t_span = (0, 20),
                              y0 = pc.vitroics,
                              method = "LSODA", t_eval = np.linspace(0, 20, 500000),# "LSODA",
                              atol = 1e-10,
                              rtol = 1e-10)


    resultsPO = sci.solve_ivp(fun = lambda t, y: po(t, y, pW = pW),
                              t_span = (0, 20),
                              y0 = pc.vitroics,
                              method = "LSODA", t_eval = np.linspace(0, 20, 100000),# "LSODA",
                              atol = 1e-10,
                              rtol = 1e-10)

    state = results.y[-1]
    #J1 = fluxes.fluxes(state, param = pc.params, ExpType = ExpType)

    print("Done S2")

    # resultsS2 = np.concatenate((np.array([resultsS2.t]), resultsS2.y)).transpose()
    # resultsS2 = pd.DataFrame(resultsS2,
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
    # resultsS2.to_csv("../results/resultsState2F1KHTest.csv")

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
    results.to_csv("../results/resultsState2KHTest.csv")

    resultsLR = np.concatenate((np.array([resultsLR.t]), resultsLR.y)).transpose()
    resultsLR = pd.DataFrame(resultsLR,
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
    resultsLR.to_csv("../results/resultsState2LRKHTest.csv")


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
    resultsPO.to_csv("../results/resultsState2POKHTest.csv")
    #print("Done State 2 calculations.")

    ## State 3 Respiration
    s2newIC = np.delete(results.tail(1).to_numpy()[0], [0])

    pc.vitroics = s2newIC
    pc.vitroics[pc.pcIS.iADP_c] = 2.5e-3

    results = sci.solve_ivp(fun = lambda t, y: s3(t, y, pW),
                            t_span = (0, 20),
                            y0 = pc.vitroics,
                            method = "Radau", t_eval = np.linspace(0, 20, 100000),#"LSODA",
                            atol = 1e-10,
                            rtol = 1e-10)


    state = results.y[-1]
    #J2 = fluxes.fluxes(state, param = pc.params, ExpType = ExpType)

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
    results.to_csv("../results/resultsState3KHTest.csv")
    #print("Done State 3 calculations.")

    ## Leak Respiration

    s2newICLR = np.delete(resultsLR.tail(1).to_numpy(), [0])

    pc.vitroics = s2newICLR
    pc.vitroics[pc.pcIS.iADP_c] = 2.5e-3

    results = sci.solve_ivp(fun = lambda t, y: lr(t, y, pW = pW),
                            t_span = (0, 20),
                            y0 = pc.vitroics,
                            method = "Radau", t_eval = np.linspace(0, 20, 500000),# "LSODA",
                            atol = 1e-10,
                            rtol = 1e-10)


    state = results.y[-1]
    #J3 = fluxes.fluxes(state, param = pc.params, ExpType = ExpType)

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
    results.to_csv("../results/resultsLeakRespKHTest.csv")
    #print("Done Leak Respiration calculations.")

    ## P/O ratio at half-maximum respiration
    s2newICPO = np.delete(resultsPO.tail(1).to_numpy()[0], [0])
    pc.vitroics = s2newICPO
    pc.vitroics[pc.pcIS.iADP_c] = 2.5e-3*0.008
    pc.vitroics[pc.pcIS.iATP_c] = 2.0e-3


    results = sci.solve_ivp(fun = lambda t, y: po(t, y, pW = pW),
                            t_span = (0, 20),
                            y0 = pc.vitroics,
                            method = "Radau", t_eval = np.linspace(0, 20, 100000),# "LSODA",
                            atol = 1e-10,
                            rtol = 1e-10)


    state = results.y[-1]
    #J4 = fluxes.fluxes(state, param = pc.params, ExpType = ExpType)

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
    results.to_csv("../results/resultsPOKHTest.csv")
    #print("Done P/O ratio calculations.")


#start = time.time()
main()
#end = time.time()
#print(end-start)

#finalConditions = np.array(a.tail(1)) ## Remove the first element before use
