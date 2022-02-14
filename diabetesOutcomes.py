import numpy as np
import pandas as pd

import pc as pc
import fluxes as fl

J_AtC = 1.256e-3 ## Used by Edwards et al.
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

## Conversion factor taken from Edwards et al 2019
rho_m = 3.0545e-6
convert = 1e9 * rho_m

def main():
    pc.params[38] = 1.15*pc.params[38]

    normal = pd.read_csv("../results/resultsDiabetes.csv").tail(1).to_numpy()[0]
    normal = np.delete(normal, [0, 1])
    flVals = fl.fluxes(normal, param=pc.params, ExpType=ExpType, potassiumW=1.0)
    JO2 = (flVals[pc.pcFl.J_C4] / 2.) * convert
    JATP = flVals[pc.pcFl.J_F1] * convert
    PO = JATP / (2. * JO2)
    ATPc = normal[pc.pcIS.iATP_c]
    ADPc = normal[pc.pcIS.iADP_c]
    ARatio = ATPc / ADPc
    pmf = normal[pc.pcIS.idPsi] - ((pc.pcPC.RT * np.log(10)) / pc.pcPC.F) * \
          ((10 ** (-normal[pc.pcIS.iH_i])) - (10 ** (-normal[pc.pcIS.iH_x])))
    values = {"JO2": JO2, "JATP": JATP, "PO": PO, "PMF": pmf, "ATPc": ATPc * 1000.,
              "ATP/ADP Ratio": ARatio}
    print(values)

main()