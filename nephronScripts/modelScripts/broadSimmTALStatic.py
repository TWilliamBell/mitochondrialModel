import itertools
import pandas as pd
import scipy.optimize as sci
import warnings as warnings

import equations
import pc

J_AtC = 4.39e-4 ## Used by Edwards et al.
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

hleaknorm = pc.params[38]
pc.pcPC.k_O2 = pc.pcPC.k_O2 / 2.0

def main():
    ics = pc.ics
    w = [0.25, 0.5, 0.75, 1.0]
    h = [1., 1.15, 5, 10]
    o2 = [0.25, 0.5, 1]
    glyc = [0.0, 1.21e-6, 9.8e-6, 1.83e-5]

    steadyStates = list()
    badRows = list()
    k = 0
    for i in itertools.product(w, w, w, w):
        for j in itertools.product(h, glyc, o2):
            print(k)
            k += 1
            pc.params[38] = list(j)[0]*hleaknorm
            ics[pc.pcIS.iO2_x] = list(j)[2]*pc.ics[pc.pcIS.iO2_x]
            f = lambda y: equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType,
                              w = list(i),
                              glyc = list(j)[1],
                              tubule = "mTAL")

            warnings.filterwarnings('error')
            try:
                row = sci.root(f, x0 = pc.finalConditions, method = "lm")
            except:
                warnings.filterwarnings('ignore')
                row = sci.root(f, x0=pc.finalConditions, method = "lm")
                badRow = list(i) + list(j)
                badRows.append(badRow)

            steadyStates.append(row)

    badRows = pd.DataFrame(badRows)
    badRows.to_csv("../results/rowsWithWarningsBroadSimmTAL.csv")
    steadyStates = pd.DataFrame(steadyStates,
                                columns= ["H_x", "dPsi", "ATP_x", "ADP_x",
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
    steadyStates.to_csv("../results/tailsBroadSimmTALStatic.csv")

main()