import numpy as np
import pandas as pd
import itertools
from collections import Counter
import pprint

import pc as pc
import equations

J_AtC = 1.234e-3#4.39e-4 ## Used by Edwards et al.
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

#pc.pcPC.k_O2 = pc.pcPC.k_O2 / 2.0

hleaknorm = pc.params[38]

def f(y): ## Differential equations, with optional arguments specified
    return equations.conservationEqs1(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType#,
                              #tubule = "mTAL"
                              )


def get_key(dict, val):
    dict = dict.__dict__
    for key in dict.keys():
        if val == dict[key]:
            return key
    return "key doesn't exist"

def main(): ## Runs differential equation for time span and outputs results to
    ## a csv file and a feather file.
    results = pd.read_csv("../results/tailsDrugSimStatic.csv")

    w = [0.25, 0.5, 0.75, 1]
    h = [1., 1.15, 1.5, 2, 5, 10]
    p = [1., 1.15, 2]

    k = 0
    l = 0
    steadyStates = list()
    badRowKeys = list()
    for i in itertools.product(w, w, w, w):
        for j in itertools.product(h, p):
            if list(j) == [1., 1.]:
                continue
            row = np.delete(results.loc[k, ].to_numpy(), [0])
            k+=1
            pc.params[38] = list(j)[0]*hleaknorm
            f = lambda y: equations.conservationEqs1(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType,
                              w = list(i),
                              DCA = list(j)[1])
            if np.linalg.norm(f(row)) < 1e-2:
                l+=1
            if np.linalg.norm(f(row)) > 1e-2:
                keys = list()
                index = np.where(np.abs(f(row)) > 1e-3)[0]
                for s in index:
                    keys.append(get_key(pc.pcIS, s))
                badRowKeys = badRowKeys + keys
    freqs = Counter(badRowKeys)
    pprint.pprint(freqs)
    print(l, flush = True)

main()