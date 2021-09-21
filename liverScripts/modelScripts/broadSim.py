import itertools
import scipy.integrate as sci
import numpy as np
import pandas as pd
import time as time

from classDefs import timeError
import equations
import pc

J_AtC = 7.0e-4 ## Used by Edwards et al.
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

pc.finalConditions[pc.pcIS.iNADH_x] = pc.pcPC.NADtot/2.0
o2norm = pc.finalConditions[pc.pcIS.iO2_x]

def main(): ## Runs differential equation for time span and outputs results to
    ## a csv file and a feather file.
    ics = pc.livFinalConditions
    w = [0.25, 0.5, 0.75, 1.0]
    o2 = [0.25, 0.5, 1]
    glyc = [0.0, 2.3e-4, 3.4e-4]

    k = 0
    for i in itertools.product(w, w, w, w):
        for j in itertools.product(glyc, o2):
            k+=1
            ics[pc.pcIS.iO2_x] = list(j)[1]*o2norm
            a = time.time()
            f = lambda t, y: equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType,
                              w = list(i),
                              timeStart = a,
                              glyc = list(j)[0])
            try:
                results = sci.solve_ivp(fun = f,
                            t_span = (0, 1000),
                            y0 = ics,
                            method = "LSODA",
                            atol = 1e-8,
                            rtol = 1e-8)
            except timeError:
                continue
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
            results.to_csv("../results/resultsBroadATP"+str(k)+".csv")
            print(results)

main()