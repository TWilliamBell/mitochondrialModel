import itertools
import scipy.optimize as sci
import pandas as pd
import warnings
import numpy as np
from collections import Counter

import equations
import pc
import classDefs

J_AtC = 1.256e-3 ## Used by Edwards et al.
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

hleaknorm = pc.params[38]

def main(): ## Runs differential equation for time span and outputs results to
    ## a csv file and a feather file.
    w = [0.25, 0.5, 0.75, 1]
    h = [1., 1.15, 1.5, 2, 5, 10]
    p = [1., 1.15, 2]

    k = 0
    steadyStates = list()
    badRows = list()
    term = list()
    for i in itertools.product(w, w, w, w):
        for j in itertools.product(h, p):
            if list(j) == [1., 1.]:
                continue
            k += 1
            print(k)
            pc.params[38] = list(j)[0]*hleaknorm
            f = lambda y: equations.conservationEqs1(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType,
                              w = list(i),
                              DCA = list(j)[1])
            try:
                result = sci.root(f, x0=pc.finalConditions, method='lm')
                row = result.x
                if result.success == False:
                    raise classDefs.didNotConverge
                else:
                    term = term + [result.message]
                    print(result.message)
            except:
                ran = np.linspace(start = 0., stop = 1., num = 10)
                guess = pc.finalConditions
                for m in ran:
                    pc.params[38] = list(j)[0] * hleaknorm * m + hleaknorm * (1 - m)
                    f = lambda y: equations.conservationEqs1(y, J_AtC=J_AtC,
                                    ExpType=ExpType,
                                    StateType=StateType,
                                    w=np.array(list(i)) * m + np.array([1, 1, 1, 1])*(1 - m),
                                    DCA=list(j)[1] * m + (1 - m))
                    results = sci.root(f, x0=guess, method='lm')
                    guess = results.x
                    term = term + [results.message]
                    print(results.message)
                row = guess
                didItWork = np.linalg.norm(f(row))
                badRow = list(i) + list(j) + list([didItWork])
                badRows.append(badRow)
                print(didItWork, flush=True)

            steadyStates.append(row)

    terms = Counter(term)
    #print(terms)

    badRows = pd.DataFrame(badRows)
    badRows.to_csv("../results/rowsWithWarningsDrugSimStatic.csv")
    steadyStates = pd.DataFrame(steadyStates,
                           columns = ["H_x", "dPsi", "ATP_x", "ADP_x",
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
    steadyStates.to_csv("../results/tailsDrugSimStatic.csv")


main()