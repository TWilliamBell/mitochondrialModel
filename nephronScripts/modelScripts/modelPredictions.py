import scipy.integrate as sci
import feather
import numpy as np
import pandas as pd

import equations
import pc

J_AtC = 1.256e-3 ## Used by Edwards et al.
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

pc.ics = pc.finalConditions

o2norm = pc.ics[pc.pcIS.iO2_x]
hleaknorm = pc.params[38]
hcnorm = pc.ics[pc.pcIS.iH_c]
kcnorm = pc.ics[pc.pcIS.iK_c]
mgcnorm = pc.ics[pc.pcIS.iMg_c]

def main(): ## Runs differential equation for time span and outputs results to
    ## a csv file and a feather file.
    for i in range(21):
        pc.ics[pc.pcIS.iO2_x] = o2norm
        pc.ics[pc.pcIS.iH_c] = hcnorm
        pc.ics[pc.pcIS.iK_c] = kcnorm
        pc.ics[pc.pcIS.iMg_c] = mgcnorm
        pc.params[38] = hleaknorm
        pW = 1.0
        J_AtC = 1.256e-3

        ## Conditions to reproduce (for the altered model) Table 3's observations from Edwards et al
        if i == 0:
            pc.ics[pc.pcIS.iO2_x] = o2norm*(60/50)
        if i == 1:
            pc.ics[pc.pcIS.iO2_x] = o2norm*(40/50)
        if i == 2:
            pc.ics[pc.pcIS.iO2_x] = o2norm*(20/50)
        if i == 3:
            pc.ics[pc.pcIS.iO2_x] = o2norm*(10/50)
        if i == 4:
            J_AtC = 0.75*J_AtC
        if i == 5:
            J_AtC = 1.25*J_AtC
        if i == 6:
            J_AtC = 1.5*J_AtC
        if i == 7:
            pc.ics[pc.pcIS.iH_c] = 10**-7.4
        if i == 8:
            pc.ics[pc.pcIS.iH_c] = 10**-7.0
        if i == 9:
            pc.ics[pc.pcIS.iH_c] = 10**-6.8
        if i == 10:
            pc.ics[pc.pcIS.iK_c] = 60.0e-3
        if i == 11:
            pc.ics[pc.pcIS.iK_c] = 140.0e-3
        if i == 12:
            pc.ics[pc.pcIS.iMg_c] = 0.2e-3
        if i == 13:
            pc.ics[pc.pcIS.iMg_c] = 0.8e-3
        if i == 14:
            pc.params[38] = 0.0
            pW = 0.0
        if i == 15:
            pc.params[38] = 0.0
        if i == 16:
            pc.params[38] = 0.0
            pW = 10.0
        if i == 17:
            pc.params[38] = hleaknorm*10.0
            pW = 0.0
        if i == 18:
            pc.params[38] = hleaknorm*10.0
        if i == 19:
            pW = 10.0
        if i == 20:
            pc.params[38] = hleaknorm*10.0
            pW = 10.0

        f = lambda t, y : equations.conservationEqs1(y, J_AtC = J_AtC, ExpType = ExpType, StateType = StateType, potassiumW = pW)

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
        results.to_csv("../results/resultsMP"+str(i)+".csv")
        feather.write_dataframe(results, "../results/resultsMP"+str(i)+"Fe.feather")
        print(results)

main()