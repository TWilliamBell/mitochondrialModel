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

def f(t, y, J_AtC = 7.0e-4, glyc = 2.3e-4,
      w = [1., 1., 1., 1.]): ## Differential equations, with optional arguments specified
    print(t)
    return equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType,
                              glyc = glyc, w = w)

## We consider cases where the PT was ATP depleted but the mTAL was not
## We see how V_mito, J_AtC, and Complex IV k_O2 impact the robustness observed
## for the mTAL.
## We aim to find necessary and sufficient conditions for the robustness
## observed in these tissues due to these three parameter changes.

def main(): ## Runs differential equation for time span and outputs results to
    ## a csv file and a feather file.

    normalO2 = pc.finalConditions[pc.pcIS.iO2_x]
    normalATP = pc.finalConditions[pc.pcIS.iATP_c]
    normalADP = pc.finalConditions[pc.pcIS.iADP_c]
    normalAMP = pc.finalConditions[pc.pcIS.iAMP_c]
    for i in range(8):
        pc.finalConditions[pc.pcIS.iO2_x] = \
            normalO2 * 0.1
        pc.finalConditions[pc.pcIS.iATP_c] = normalATP
        pc.finalConditions[pc.pcIS.iADP_c] = normalADP
        pc.finalConditions[pc.pcIS.iAMP_c] = normalAMP
        if i == 0: ## Finished this on my computer
            g = lambda t, y: f(t, y, J_AtC = 1.256e-3*1.25) ## J_AtC
        if i == 1:
            g = lambda t, y: f(t, y, glyc = 0.0) ## glycolysis
        if i == 2:
            g = lambda t, y: f(t, y, J_AtC = 1.256e-3*1.25, glyc = 0.0)
            ## J_AtC, glycolysis
        if i == 3:
            pc.finalConditions[pc.pcIS.iATP_c] = 0.0021
            pc.finalConditions[pc.pcIS.iADP_c] = 0.000527
            pc.finalConditions[pc.pcIS.iAMP_c] = 8.573543e-05
            g = lambda t, y: f(t, y, J_AtC = 1.256e-3*1.25, glyc = 0.0)
        if i == 4:
            g = lambda t, y: f(t, y, J_AtC = 1.256e-3*1.25, glyc = 0.0,
                               w = [0.75/1.2, 0.5/0.525, 0.25/0.4, 1.])
        if i == 5:
            g = lambda t, y: f(t, y,
                               w = [0.75/1.2, 0.5/0.525,
                                    0.25/0.4, 1.]) ## Just OXPHOS
        if i == 6:
            g = lambda t, y: f(t, y, J_AtC = 1.256e-3*1.25,
                               w=[0.75 / 1.2, 0.5 / 0.525,
                                  0.25 / 0.4, 1.])
        if i == 7:
            g = lambda t, y: f(t, y, glyc = 0.0,
                               w=[0.75 / 1.2, 0.5 / 0.525,
                                  0.25 / 0.4, 1.])


        results = sci.solve_ivp(fun = g,
                            t_span = (0, 120),
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
        results.to_csv("../results/resultsMech"+str(i+1)+".csv")
        #feather.write_dataframe(results, "../results/resultsMech"+str(i+1)+
        #                        ".feather")
        print(results)

    pc.finalConditions[pc.pcIS.iO2_x] = normalO2
    pc.finalConditions[pc.pcIS.iATP_c] = normalATP
    pc.finalConditions[pc.pcIS.iADP_c] = normalADP
    pc.finalConditions[pc.pcIS.iAMP_c] = normalAMP
    for i in range(8):
        pc.finalConditions[pc.pcIS.iATP_c] = normalATP
        pc.finalConditions[pc.pcIS.iADP_c] = normalADP
        pc.finalConditions[pc.pcIS.iAMP_c] = normalAMP
        if i == 0: ## Finished on cpu150
            g = lambda t, y: f(t, y, J_AtC = 1.256e-3,
                               w = [1., 0.25, 1., 1.]) ## J_AtC
        if i == 1:
            g = lambda t, y: f(t, y, glyc = 0.0,
                               w = [1., 0.25, 1., 1.]) ## glycolysis
        if i == 2:
            g = lambda t, y: f(t, y, J_AtC = 1.256e-3, glyc = 0.0,
                               w = [1., 0.25, 1., 1.])
            ## J_AtC, glycolysis
        if i == 3:
            pc.finalConditions[pc.pcIS.iATP_c] = 0.0021
            pc.finalConditions[pc.pcIS.iADP_c] = 0.000527
            pc.finalConditions[pc.pcIS.iAMP_c] = 8.573543e-05
            g = lambda t, y: f(t, y, J_AtC = 1.256e-3, glyc = 0.0,
                               w = [1., 0.25, 1., 1.])
        if i == 4:
            g = lambda t, y: f(t, y, J_AtC = 1.256e-3, glyc = 0.0,
                               w = [0.75/1.2, 0.25*0.5/0.525,
                                    0.25/0.4, 1.]) ## Both+OXPHOS
        if i == 5:
            g = lambda t, y: f(t, y,
                               w = [0.75/1.2, 0.25*0.5/0.525,
                                    0.25/0.4, 1.]) ## Just OXPHOS
        if i == 6:
            g = lambda t, y: f(t, y, J_AtC=1.256e-3,
                                w = [0.75 / 1.2, 0.25 * 0.5 / 0.525,
                                        0.25 / 0.4, 1.])
        if i == 7:
            g = lambda t, y: f(t, y, glyc=0.0,
                                w = [0.75 / 1.2, 0.25 * 0.5 / 0.525,
                                    0.25 / 0.4, 1.])


        results = sci.solve_ivp(fun = g,
                            t_span = (0, 120),
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
        results.to_csv("../results/resultsMech"+str(i+1)+"C3.csv")
        print(results)

main()
