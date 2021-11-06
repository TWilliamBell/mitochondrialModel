import itertools
import scipy.integrate as sci
import numpy as np
import pandas as pd
import feather

import fluxes as fl
import pc

J_AtC = 1.256e-3 ## Used by Edwards et al.
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

hleaknorm = pc.params[38]


def main():
    drugSimDat = pd.read_csv(
        "../results/tailsDrugSim.csv")
    drugSimDat = drugSimDat.sort_values(by = "nums",
                                        ignore_index = True)
    iterProd = pd.read_feather(
        "../results/iterProdDrugSim.feather")
    iterProdMitDis = pd.read_feather(
        "../results/iterProd.feather")
    mitDisDat = pd.read_csv(
        "../results/tailsMitDisPT.csv")
    mitDisDat = mitDisDat.sort_values(by = "nums",
                                      ignore_index = True)


    fs = list()
    conditions = list()
    for i in range(len(drugSimDat)):
        y = np.delete(drugSimDat.loc[i, ].to_numpy(),
                          [0, 1, 2, 3])
        k = drugSimDat["nums"][i]
        w = iterProd.loc[k-1, ][[1, 2, 3, 4]].to_numpy()
        j = iterProd.loc[k-1, ][[5, 6]].to_numpy()
        if j[0] == 1. and j[1] == 1.:
            continue
        conditions.append(list(w)+list(j))
        pc.params[38] = list(j)[0]*hleaknorm
        f = fl.fluxes(y,
                      ExpType = ExpType,
                      param = pc.params,
                      w = w,
                      DCA = list(j)[1])
        fs.append(list(f))

    pc.params[38] = hleaknorm
    for i in range(len(mitDisDat)):
        y = np.delete(mitDisDat.loc[i, ].to_numpy(),
                          [0, 1, 2, 3])
        k = mitDisDat["nums"][i]
        w = iterProdMitDis.loc[k-1, ].to_numpy()
        conditions.append(list(w)+[1., 1.])
        f = fl.fluxes(y,
                      ExpType = ExpType,
                      param = pc.params,
                      w = w)
        fs.append(list(f))

    rho_m = 3.0545e-6
    convert = 1e9 * rho_m
    po = list()
    qo = list()
    ga = list()
    q = 0
    for i in range(len(fs)):
        ## Do some analysis here
        info = fs[i]
        ## P/O ratio
        po.append(info[pc.pcFl.J_F1]/info[pc.pcFl.J_C4])
        ## Q_O2
        qo.append(convert*info[pc.pcFl.J_C4]/2.)
        ## ATP Gen
        ga.append(convert*info[pc.pcFl.J_F1])
        # if ga[i] < 2.0:
        #     print(conditions[i])
        #     print(ga[i])
        if qo[i] > 4.88 and conditions[i][5] == 1.0:
            #print(conditions[i])
            #print(qo[i])
            q+=1
        ## Report univar info
        # if np.sum((1.-np.array(conditions[i])) == 0.) == 5:
        #     print(conditions[i])
        #     print(po[-1])
        #     print(qo[-1])
        #     print(ga[-1])
        #     print("\n")

    print("PO Ratio range:")
    print((round(min(po), 2), round(max(po), 2)))
    print("Q_O2 range:")
    print((round(min(qo), 2), round(max(qo), 2)))
    print("ATP generation range:")
    print((round(min(ga), 2), round(max(ga), 2)))
    print(q)
    #import pprint as pp
    #pp.pprint(list(np.round(qo, 2)))
    #print(np.sum(np.round(ga, 2) < 5.03)/len(ga))

main()