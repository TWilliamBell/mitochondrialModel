import feather
import numpy as np
import pandas as pd

import fluxes as flx
import pc

ExpType = 1

def main():

    dfList = list()
    ## Collect the activities of Complex I, Complex III, Complex IV,
    ## and F0F1-ATPase
    for i in range(5):
        for j in range(5):
            if i == 2 & j==2:
                dat = pd.read_feather("../results/resultsATP.feather")\
                    .tail(1).values[0]
                dat = np.delete(dat, 0)
            else:
                dat = pd.read_feather("../results/resultsIOATP"+str(i)+str(j)+
                                           "Fe.feather").tail(1).values[0]
                dat = np.delete(dat, 0)
            dat = flx.fluxes(dat, pc.params, ExpType)
            dfList.append(dat[[pc.pcFl.J_C1, pc.pcFl.J_C3,
                               pc.pcFl.J_C4, pc.pcFl.J_F1]]/dat[pc.pcFl.J_cits])
    results = pd.DataFrame(np.array(dfList), columns = ["J_C1", "J_C3", "J_C4",
                                                        "F0F1-ATPase"])
    print(results)
    print("\n")
    modelPrediction = dfList[8] ## The eighth is the baseline
    print("The ratio predicted by the baseline model:")
    print(modelPrediction)
    feather.write_dataframe(results, "../results/ratiosHealthyATP.feather")

main()