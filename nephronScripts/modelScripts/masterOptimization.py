import scipy.integrate as sci
import scipy.optimize as sco
import numpy as np
import pandas as pd

import equations
import pc
import fluxes
import fluxesmTAL

StateType = 1

# PT
def kh(t, y, pW):
    print(t)
    return equations.conservationEqs(y, J_AtC = 1.256e-3,
                              ExpType = 1,
                              StateType = StateType,
                              w = [1., 1., 1., 1.], potassiumW = pW)

def s2(t, y, pW): ## For State 2 Resp.
    print(t)
    a = equations.conservationEqs(y, J_AtC = 0.,
                              ExpType = 2,
                              StateType = StateType,
                              w = [1., 1., 1., 0.], potassiumW = pW)
    return a

def s3(t, y, pW): ## Differential equations, with optional arguments specified
    print(t)
    return equations.conservationEqs(y, J_AtC = 0.,
                              ExpType = 2,
                              StateType = StateType, potassiumW = pW)

def lr(t, y, pW): ## Differential equations, with optional arguments specified
    print(t)
    a = equations.conservationEqs(y, J_AtC = 0.,
                              ExpType = 2,
                              StateType = StateType,
                              w = [1., 1., 1., 0.], potassiumW = pW)
    a[pc.pcIS.iADP_c] = 0
    a[pc.pcIS.iATP_c] = 0
    return a

def po(t, y, pW):
    print(t)
    a = equations.conservationEqs(y, J_AtC = 0.,
                              ExpType = 2,
                              StateType = StateType,
                              w = [1., 1., 1., 1.], potassiumW = pW)
    a[pc.pcIS.iADP_c] = 0
    a[pc.pcIS.iATP_c] = 0
    return a


## mTAL
def s2mTAL(t, y, pW): ## For State 2 Resp.
    print(t)
    a = equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = 2,
                              StateType = StateType,
                              w = [1., 1., 1., 0.],
                              tubule = "mTAL", potassiumW = pW)
    return a

def s3mTAL(t, y, pW): ## Differential equations, with optional arguments specified
    print(t)
    return equations.conservationEqs(y, J_AtC = 0.,
                              ExpType = 2,
                              StateType = StateType,
                              tubule = "mTAL", potassiumW = pW)

def pomTAL(t, y, pW):
    print(t)
    a = equations.conservationEqs(y, J_AtC = 0.,
                              ExpType = 2,
                              StateType = StateType,
                              tubule = "mTAL", potassiumW = pW)
    a[pc.pcIS.iADP_c] = 0
    a[pc.pcIS.iATP_c] = 0
    return a

def cost(PTpred, mTALpred, ko2pt, ko2mTAL):
    S3 = PTpred[0]
    RCR = PTpred[1]
    LR = PTpred[2]
    PO = PTpred[3]
    RCRmTAL = mTALpred[0]
    POmTAL = mTALpred[1]
    if PO > POmTAL:
        compare = 1
    else:
        compare = 0
    if ko2pt < ko2mTAL:
        compare+=1

    return ((S3-4.08)/0.6)**2+((RCR-8.5)/0.9)**2+((LR-0.37)/0.06)**2+\
           ((PO-1.8)/0.1)**2+((KH-1.05)/0.025)**2+((RCRmTAL-10.0)/1.6)**2+\
           ((POmTAL-1.92)/0.13)**2+compare

def optimFn(optParam): ## Runs differential equation for time span and outputs results to
    ## a csv file and a feather file.
    khpt = optParam[0] ## Done
    khmTAL = optParam[1]
    ko2pt = optParam[2] ## Done - need to test global
    ko2mTAL = optParam[3]
    hleakpt = optParam[4] ## Not done
    hleakmTAL = optParam[5]
    kleakpt = optParam[6] ## Done - need to test args
    kleakmTAL = optParam[7]
    antpt = optParam[8]
    antmTAL = optParam[9]

    pW = kleakpt

    ## Proximal Tubule Predictions
    global pc
    pc.pcPC.k_O2 = 1.2e-4*ko2pt
    ## The in vivo case
    ## Increase Potassium-Hydrogen Antiporter activity by 100x from 'new normal'
    ## and see if we get that consistent 1.05 increase in dPsi

    pc.params[37] = khpt * (4.7580e+06 / 15) * 20
    pc.params[39] = 347.4 * hleakpt
    pc.params[34] = 0.00675 * antpt

    resultsKH = sci.solve_ivp(fun = lambda t, y: kh(t, y, pW = pW),
                              t_span = (0, 10),
                              y0 = pc.ics,
                              method = "LSODA",
                              atol = 1e-10,
                              rtol = 1e-10)

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
    pc.params[37] = khpt * 4.7580e+06 / 15

    resultsKHbase = sci.solve_ivp(fun = lambda t, y: kh(t, y, pW = pW),
                              t_span = (0, 10),
                              y0 = pc.ics,
                              method = "LSODA",
                              atol = 1e-10,
                              rtol = 1e-10)

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

    pc.params[37] = khpt * 4.7580e+06 / 15

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

    resultsS2 = sci.solve_ivp(fun = lambda t, y: s2(t, y, pW = pW),
                              t_span = (0, 20),
                              y0 = pc.vitroics,
                              method = "LSODA",
                              atol = 1e-10,
                              rtol = 1e-10)

    results = sci.solve_ivp(fun = lambda t, y: s3(t, y, pW = pW),
                            t_span = (0, 20),
                            y0 = pc.vitroics,
                            method = "LSODA",
                            atol = 1e-10,
                            rtol = 1e-10)

    resultsLR = sci.solve_ivp(fun = lambda t, y: lr(t, y, pW = pW),
                              t_span = (0, 20),
                              y0 = pc.vitroics,
                              method = "LSODA",
                              atol = 1e-10,
                              rtol = 1e-10)


    resultsPO = sci.solve_ivp(fun = lambda t, y: po(t, y, pW = pW),
                              t_span = (0, 20),
                              y0 = pc.vitroics,
                              method = "LSODA",
                              atol = 1e-10,
                              rtol = 1e-10)

    state = results.y[-1]
    #J1 = fluxes.fluxes(state, param = pc.params, ExpType = ExpType)


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
    resultsS2.to_csv("../results/resultsState2F1KHTest.csv")

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
    print("Done State 2 calculations.")

    ## State 3 Respiration
    s2newIC = np.delete(results.tail(1).to_numpy()[0], [0])

    pc.vitroics = s2newIC
    pc.vitroics[pc.pcIS.iADP_c] = 2.5e-3

    results = sci.solve_ivp(fun = lambda t, y: s3(t, y, pW = pW),
                            t_span = (0, 20),
                            y0 = pc.vitroics,
                            method = "LSODA",
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
    print("Done State 3 calculations.")

    ## Leak Respiration

    s2newICLR = np.delete(resultsLR.tail(1).to_numpy(), [0])

    pc.vitroics = s2newICLR
    pc.vitroics[pc.pcIS.iADP_c] = 2.5e-3

    results = sci.solve_ivp(fun = lambda t, y: lr(t, y, pW = pW),
                            t_span = (0, 20),
                            y0 = pc.vitroics,
                            method = "LSODA",
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
    print("Done Leak Respiration calculations.")

    ## P/O ratio at half-maximum respiration
    s2newICPO = np.delete(resultsPO.tail(1).to_numpy()[0], [0])
    pc.vitroics = s2newICPO
    pc.vitroics[pc.pcIS.iADP_c] = 2.5e-3*0.008
    pc.vitroics[pc.pcIS.iATP_c] = 2.0e-3


    results = sci.solve_ivp(fun = lambda t, y: po(t, y, pW = pW),
                            t_span = (0, 20),
                            y0 = pc.vitroics,
                            method = "LSODA",
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
    print("Done P/O ratio calculations.")

    ## Calculate predictions

    rho_m = 3.0545e-6
    convert = 1e9 * rho_m

    khUp = pd.read_csv("../results/resultsKH100Test.csv").tail(1).to_numpy()[0]
    khUp = np.delete(khUp, [0, 1])

    reg = pd.read_csv("../results/resultsKHbaseTest.csv").tail(1).to_numpy()[0]
    reg = np.delete(reg, [0, 1])

    ## Data loading and cleaning, based on measureVars.py simulations
    s2dat = pd.read_csv("../results/resultsState2KHTest.csv").tail(1).to_numpy()[0]
    s2dat= np.delete(s2dat, [0, 1])
    #s2LR = pd.read_csv("../results/resultsState2LR.csv").tail(1).to_numpy()[0]
    #s2LR = np.delete(s2LR, [0, 1]


    # print("State 2 State Variables:")
    # print(s2)
    # print("\n")

    s3dat = pd.read_csv("../results/resultsState3KHTest.csv")

    # print("State 3 State Variables:")
    # print(np.delete(s3.tail(1).to_numpy()[0], [0, 1]))
    # print("\n")

    s3dat = s3dat.to_numpy()

    lrdat = pd.read_csv("../results/resultsLeakRespKHTest.csv")

    # print("Leak respiration State Variables:")
    # print(np.delete(lr.tail(1).to_numpy()[0], [0, 1]))
    # print("\n")

    lrdat = lrdat.to_numpy()

    podat = pd.read_csv("../results/resultsPOKHTest.csv")

    # print("P/O measurement State Variables:")
    # print(np.delete(po.tail(1).to_numpy()[0], [0, 1]))
    # print("\n")

    podat = podat.to_numpy()

    ## Calculate fluxes
    #Js2 = list()
    Js3 = list()
    Jlr = list()
    Jpo = list()

    #for i in range(len(s2)):
    #    s2t = np.delete(s2[i], [0, 1])
    #    Js2.append(fluxes.fluxes(s2t, param = pc.params, ExpType = ExpType))

    for i in range(len(s3dat)):
        s3t = np.delete(s3dat[i], [0, 1])
        Js3.append(fluxes.fluxes(s3t, param = pc.params, ExpType = 2,
                              w = [1., 1., 1., 1.]))

    for i in range(len(lrdat)):
        lrt = np.delete(lrdat[i], [0, 1])
        Jlr.append(fluxes.fluxes(lrt, param = pc.params, ExpType = 2,
                              w = [1., 1., 1., 1.]))

    for i in range(len(podat)):
        pot = np.delete(podat[i], [0, 1])
        Jpo.append(fluxes.fluxes(pot, param = pc.params, ExpType = 2,
                              w = [1., 1., 1., 1.]))

    #s2JO2 = [item[2] for item in Js2]
    #JO2s2 = max(s2JO2)
    #ind2 = np.argmax(s2JO2)
    #s2ATP = [item[3] for item in Js2]
    #JATPs2 = s2ATP[ind2]
    Js2 = fluxes.fluxes(s2dat, param = pc.params, ExpType = 2,
                              w = [1., 1., 1., 0.])
    JO2s2 = Js2[2]
    JATPs2 = Js2[3]

    #Js2LR = fluxes.fluxes(s2LR, param = pc.params, ExpType = ExpType)
    #JO2s2LR = Js2LR[2]
    #print("Leak respiration state 2 oxygen consumption:")
    #print(JO2s2LR)

    s3JO2 = [item[2] for item in Js3]
    JO2s3 = max(s3JO2)
    ind3 = np.argmax(s3JO2)
    s3ATP = [item[3] for item in Js3]
    JATPs3 = s3ATP[ind3]


    lrJO2 = [item[2] for item in Jlr]
    JO2lr = max(lrJO2)
    indlr = np.argmax(lrJO2)
    lrATP = [item[3] for item in Jlr]
    JATPlr = lrATP[indlr]

    poJO2 = [item[2] for item in Jpo]
    JO2po = max(poJO2)
    indpo = np.argmax(poJO2)
    poATP = [item[3] for item in Jpo]
    JATPpo = poATP[indpo]



    ## Extract values of interest
    valsdPsi = {"Up" : khUp[1], "Base" : reg[1]}
    valsS2 = {"JO2" : (JO2s2/2.)*convert, "JATP" : JATPs2*convert}
    valsS3 = {"JO2" : (JO2s3/2.)*convert, "JATP" : JATPs3*convert}
    valsLR = {"JO2" : (JO2lr/2.)*convert, "JATP" : JATPlr*convert}
    valsPO = {"JO2" : (JO2po/2.)*convert, "JATP" : JATPpo*convert}

    ## Find RC and PO ratios
    PO = valsPO["JATP"]/(2*valsPO["JO2"])
    RCR = valsS3["JO2"]/valsS2["JO2"]

    ## Difference between "Nigericin" and baseline for dPsi
    KH = valsdPsi["Up"]/valsdPsi["Base"]

    ## All values used in fitting
    PTpred = (KH, valsS3["JO2"], RCR, valsLR["JO2"], PO)

    ## mTAL predictions
    pc.pcPC.k_O2 = 1.2e-4 * ko2mTAL
    pc.params[37] = khmTAL * 4.7580e+06 / 15
    pc.params[39] = 347.4 * hleakmTAL
    pc.params[34] = 0.00675 * antmTAL

    pW = kleakmTAL

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

    results = sci.solve_ivp(fun = lambda t, y: s3mTAL(t, y, pW = pW),
                            t_span = (0, 20),
                            y0 = pc.vitroics,
                            method = "LSODA",
                            atol = 1e-10,
                            rtol = 1e-10)

    resultsPO = sci.solve_ivp(fun = lambda t, y: pomTAL(t, y, pW = pW),
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
    print("Done State 2 calculations.")

    ## State 3 Respiration
    s2newIC = np.delete(results.tail(1).to_numpy()[0], [0])

    pc.vitroics = s2newIC
    pc.vitroics[pc.pcIS.iADP_c] = 2.5e-3

    results = sci.solve_ivp(fun = lambda t, y: s3mTAL(t, y, pW = pW),
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
    print("Done State 3 calculations.")

    ## P/O ratio at half-maximum respiration
    s2newICPO = np.delete(resultsPO.tail(1).to_numpy()[0], [0])
    pc.vitroics = s2newICPO
    pc.vitroics[pc.pcIS.iADP_c] = 2.5e-3*0.008
    pc.vitroics[pc.pcIS.iATP_c] = 2.0e-3


    results = sci.solve_ivp(fun = lambda t, y: pomTAL(t, y, pW = pW),
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
    print("Done P/O ratio calculations.")

    rho_m = 3.0545e-6
    convert = 1e9 * rho_m

    ## Data loading and cleaning, based on measureVars.py simulations
    s2dat = pd.read_csv("../results/resultsState2TAL.csv").tail(1).to_numpy()[0]
    s2dat = np.delete(s2dat, [0, 1])

    s3dat = pd.read_csv("../results/resultsState3TAL.csv")
    s3dat = s3dat.to_numpy()

    podat = pd.read_csv("../results/resultsPOTAL.csv")
    podat = podat.to_numpy()

    ## Calculate fluxes
    Js3 = list()
    Jpo = list()

    for i in range(len(s3dat)):
        s3t = np.delete(s3dat[i], [0, 1])
        Js3.append(fluxes.fluxesmTAL(s3t, param = pc.params, ExpType = 2))

    for i in range(len(podat)):
       pot = np.delete(podat[i], [0, 1])
       Jpo.append(fluxes.fluxesmTAL(pot, param = pc.params, ExpType = 2))

    Js2 = fluxes.fluxesmTAL(s2dat, param = pc.params, ExpType = 2)
    JO2s2 = Js2[2]
    JATPs2 = Js2[3]

    s3JO2 = [item[2] for item in Js3]
    JO2s3 = max(s3JO2)
    ind3 = np.argmax(s3JO2)
    s3ATP = [item[3] for item in Js3]
    JATPs3 = s3ATP[ind3]

    poJO2 = [item[2] for item in Jpo]
    JO2po = max(poJO2)
    indpo = np.argmax(poJO2)
    poATP = [item[3] for item in Jpo]
    JATPpo = poATP[indpo]

    ## Extract values of interest
    valsS2 = {"JO2" : (JO2s2/2.)*convert, "JATP" : JATPs2*convert}
    valsS3 = {"JO2" : (JO2s3/2.)*convert, "JATP" : JATPs3*convert}
    valsPO = {"JO2" : (JO2po/2.)*convert, "JATP" : JATPpo*convert}

    ## Find RC and PO ratios
    PO = valsPO["JATP"]/(2.*valsPO["JO2"])
    RCR = valsS3["JO2"]/valsS2["JO2"]

    ## All values used in fitting
    mTALpred = (valsS2["JO2"], valsS3["JO2"], RCR, PO)

    ## Compute cost
    cost = cost(PTpred, mTALpred, ko2pt, ko2mTAL)

def main():

    results = sco.minimize(fun = optimFn,
                           x0 = np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]),
                           bounds = ((0.1, 2.),
                                     (0.1, 2.),
                                     (0.1, 2.),
                                     (0.1, 2.),
                                     (0.1, 2.),
                                     (0.1, 2.),
                                     (0.1, 2.),
                                     (0.1, 2.),
                                     (0.1, 2.),
                                     (0.1, 2.)))
    a = pd.DataFrame(results.x)
    a.to_csv("../results/optimizationOutput.csv")
    print(a.success)

main()