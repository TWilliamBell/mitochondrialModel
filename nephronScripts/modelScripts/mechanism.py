import scipy.integrate as sci
import feather
import numpy as np
import pandas as pd
#import time

import equations
import pc

ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

def f(t, y, J_AtC = 1.7e-3, w = [1., 1., 1., 1.]): ## Differential equations, with optional arguments specified
    return equations.conservationEqsmTAL(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType,
                              w = w)

def h(t, y, J_AtC = 1.7e-3, w = [1., 1., 1., 1.]):
    return equations.conservationEqs1(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType, params = pc.paramsmTAL,
                              tubule = "PT", w = w)

## We consider cases where the PT was ATP depleted but the mTAL was not
## We see how V_mito, J_AtC, and Complex IV k_O2 impact the robustness observed
## for the mTAL.
## We aim to find necessary and sufficient conditions for the robustness
## observed in these tissues due to these three parameter changes.

def main(): ## Runs differential equation for time span and outputs results to
    ## a csv file and a feather file.
    normalO2 = pc.finalConditions[pc.pcIS.iO2_x]
    normalKO2 = pc.pcPC.k_O2

    for i in range(6):
        pc.pcPC.k_O2 = normalKO2 / 2.0
        pc.finalConditions[pc.pcIS.iO2_x] = \
            normalO2 * 0.1
        if i == 0:
            pc.pcPC.k_O2 = normalKO2 ## k_O2
            g = lambda t, y: f(t, y)
        if i == 1:
            g = lambda t, y: f(t, y, J_AtC = 1.256e-3) ## J_AtC
        if i == 2:
            pc.finalConditions[pc.pcIS.iO2_x] = 0.1*(normalO2 / 5.0)
            g = lambda t, y: h(t, y) ## V_mito
        if i == 3:
            pc.finalConditions[pc.pcIS.iO2_x] = 0.1*(normalO2 / 5.0)
            g = lambda t, y: h(t, y, J_AtC = 1.256e-3)
            ## V_mito, J_AtC
        if i == 4:
            pc.pcPC.k_O2 = normalKO2
            g = lambda t, y: f(t, y, J_AtC = 1.256e-3)
            ## J_AtC, k_O2
        if i == 5:
            pc.pcPC.k_O2 = normalKO2
            pc.finalConditions[pc.pcIS.iO2_x] = 0.1*(normalO2 / 5.0)
            g = lambda t, y: h(t, y) ## V_mito, k_O2

        results = sci.solve_ivp(fun = g,
                            t_span = (0, 100000),
                            y0 = pc.finalConditions,
                            method = "LSODA",
                            atol = 1e-9,
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
        results.to_csv("../results/resultsmTALMech"+str(i+1)+".csv")
        print(results)

    pc.finalConditions[pc.pcIS.iO2_x] = normalO2
    for i in range(4):
        pc.pcPC.k_O2 = normalKO2 / 2.0
        pc.finalConditions[pc.pcIS.iO2_x] = normalO2
        if i == 0:
            pc.pcPC.k_O2 = normalKO2 ## k_O2
            g = lambda t, y: f(t, y, w = [1., 0.25, 1., 1.])
        if i == 1:
            g = lambda t, y: f(t, y, J_AtC = 1.256e-3,
                               w = [1., 0.25, 1., 1.]) ## J_AtC
        if i == 2:
            pc.finalConditions[pc.pcIS.iO2_x] = (normalO2 / 5.0)
            g = lambda t, y: h(t, y, w = [1., 0.25, 1., 1.]) ## V_mito
        if i == 3:
            pc.finalConditions[pc.pcIS.iO2_x] = (normalO2 / 5.0)
            g = lambda t, y: h(t, y, J_AtC = 1.256e-3,
                               w = [1., 0.25, 1., 1.])
            ## V_mito, J_AtC
        # if i == 4:
        #     pc.pcPC.k_O2 = normalKO2
        #     pc.finalConditions[pc.pcIS.iO2_x] = (normalO2 / 5.0)
        #     g = lambda t, y: h(t, y, J_AtC=1.256e-3,
        #                        w=[1., 0.25, 1., 1.])
        #     ## V_mito, J_AtC, k_O2

        results = sci.solve_ivp(fun = g,
                            t_span = (0, 100000),
                            y0 = pc.finalConditions,
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
        results.to_csv("../results/resultsmTALMech"+str(i+1)+"C3.csv")
        print(results)

#start = time.time()
main()
#end = time.time()
#print(end-start)

#finalConditions = np.array(a.tail(1)) ## Remove the first element before use

