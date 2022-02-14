import numpy as np
import numdifftools as nd
import sympy as sym

import equations
from classDefs import outOfOrder
import pc

print(np.array(pc.finalConditions)-np.array(pc.ics))

J_AtC = 1.256e-3 ## Used by Edwards et al.
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

def dL():
    return np.array([pc.pcIS.iG6P_c, ## Confirmed directly
                     pc.pcIS.iGLC_c, ## Confirmed directly
                     pc.pcIS.iFUM_c, ## Confirmed directly
                     pc.pcIS.iFUM_i, ## Confirmed directly
                     pc.pcIS.iK_i, ## Confirmed directly
                     pc.pcIS.iMg_i, ## Confirmed directly
                     pc.pcIS.iH_i ## Confirmed directly
                     ], dtype = int)

def chooseK(i, k, delList = dL()):
    if not all(delList[i] >= delList[i+1] for i in range(len(delList)-1)):
        raise outOfOrder() ## Performs check for what follows that uses delList
    for j in range(len(delList)):
        if k + 1 >= delList[j]:
            k = i + len(delList) - j
            break
    if i < np.min(delList):
        k = i
    return k

def f(y): ## Differential equations, with optional arguments specified, and many states as well
    x = pc.finalConditions
    k = -1
    for i in range(len(y)): ## So that we do not have to provide every state
        ## Depends heavily on the indices being in decreasing order
        k = chooseK(i, k)
        if y[i] != x[k]: ## If the state provided differs from the default final condition,
            ## prefer the value provided
            x[k] = y[i]
    result = equations.conservationEqs(x, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType)
    result = np.delete(result, dL()) ## For the linearization, we only wish to provide
    ## the states that we've not excluded
    return result

def main():
    jac = nd.Jacobian(f)
    equi = np.array(pc.finalConditions)
    ## List of state variables we delete from final condition calculation
    equi = np.delete(equi, dL())
    print(dL())
    vals, vecs = np.linalg.eig(jac(equi))
    print(vals)

main()