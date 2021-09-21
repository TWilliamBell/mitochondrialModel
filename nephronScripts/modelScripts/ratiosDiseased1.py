import feather
import numpy as np
import pandas as pd
import itertools

import fluxes as flx
import pc

ExpType = 1

def main():

    dfList = list()
    ## Collect the activities of Complex I, Complex III, Complex IV,
    ## and F0F1-ATPase
    w = [0.25, 0.5, 0.75, 1]
    i = 0
    for j in itertools.product(w, w, w, w):
        i+=1
        try:
            dat = pd.read_feather("../results/resultsDiseaseATP"+str(i+1)+
                                       "Fe.feather").tail(1).values[0]
        except:
            continue
        dat = np.delete(dat, 0)
        dat = flx.fluxes(dat, pc.params, ExpType, w = list(j))
        dfList.append(dat[[pc.pcFl.J_C1, pc.pcFl.J_C3,
                           pc.pcFl.J_C4, pc.pcFl.J_F1]]/dat[pc.pcFl.J_cits])
    results = pd.DataFrame(np.array(dfList), columns = ["J_C1", "J_C3", "J_C4",
                                                        "F0F1-ATPase"])
    print(results)
    print("\n")
    modelPrediction = dfList[252] ## The eighth is the baseline
    print("The ratio predicted by the baseline model:")
    print(modelPrediction)
    feather.write_dataframe(results, "../results/ratiosDiseasedATP.feather")

main()