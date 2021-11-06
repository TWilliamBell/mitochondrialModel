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

def f(t, y, i = 1., w = [1., 1., 1., 1.]): ## Differential equations, with optional arguments specified
    return equations.conservationEqs1(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType,
                              capacitanceWeight = 1./i,
                              w = w) ## Choose how much I want
                                ## to reduce the capacitance.

def main(): ## Runs differential equation for time span and outputs results to
    ## a csv file and a feather file.
    for i in range(2, 10):
        ## The capacitance changed with healthy cells
        # g = lambda t, y: f(t, y, i = i) ## Choosing a range of values for
        # ## capacitance
        # results = sci.solve_ivp(fun = g,
        #                     t_span = (0, 100000),
        #                     y0 = pc.ics,
        #                     method = "LSODA",
        #                     atol = 1e-8,
        #                     rtol = 1e-10)
        ## The capacitance changed with cyanide cells
        h = lambda t, y: f(t, y, i = i, w = [1., 1., 0.00048, 1.])
        results2 = sci.solve_ivp(fun = h,
                            t_span = (0, 100000),
                            y0 = pc.ics,
                            method = "LSODA",
                            atol = 1e-8,
                            rtol = 1e-10)
        # results = np.concatenate((np.array([results.t]), results.y)).transpose()
        # results = pd.DataFrame(results,
        #                    columns = ["t", "H_x", "dPsi", "ATP_x", "ADP_x",
        #                               "AMP_x", "GTP_x", "GDP_x", "Pi_x", "NADH_x",
        #                               "QH2_x", "OAA_x", "ACCOA_x", "CIT_x", "ICIT_x",
        #                               "AKG_x", "SCOA_x", "COASH_x", "SUC_x", "FUM_x",
        #                               "MAL_x", "GLU_x", "ASP_x", "K_x", "Mg_x",
        #                               "O2_x", "CO2tot_x", "Cred_i", "ATP_i", "ADP_i",
        #                               "AMP_i", "Pi_i", "H_i", "Mg_i", "K_i", "ATP_c",
        #                               "ADP_c", "Pi_c", "H_c", "Mg_c", "K_c", "PYR_x",
        #                               "PYR_i", "PYR_c", "CIT_i", "CIT_c", "AKG_i",
        #                               "AKG_c", "SUC_i", "SUC_c", "MAL_i", "MAL_c",
        #                               "ASP_i", "ASP_c", "GLU_i", "GLU_c", "FUM_i",
        #                               "FUM_c", "ICIT_i", "ICIT_c", "GLC_c", "G6P_c",
        #                               "PCr_c", "AMP_c"])
        # results.to_csv("../results/resultsCap"+str(i)+".csv")
        # feather.write_dataframe(results, "../results/resultsCap"+str(i)+".feather")

        results2 = np.concatenate((np.array([results2.t]), results2.y)).transpose()
        results2 = pd.DataFrame(results2,
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
        results2.to_csv("../results/resultsCapCyan"+str(i)+".csv")
        feather.write_dataframe(results2, "../results/resultsCapCyan"+str(i)+".feather")

#start = time.time()
a = main()
#end = time.time()
#print(end-start)

print(a)
#finalConditions = np.array(a.tail(1)) ## Remove the first element before use