import numpy as np
import numdifftools as nd
import pseudopy as ps
from matplotlib import pyplot

import equations
from classDefs import outOfOrder
import pc

J_AtC = 1.256e-3 ## Used by Edwards et al.
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

def chooseK(i, k, delList):
    for j in range(len(delList)):
        if k + 1 >= delList[j]:
            k = i + len(delList) - j
            break
    if i < np.min(delList):
        k = i
    return k

def f(y): ## Differential equations, with optional arguments specified,
    ## and many states as well
    x = pc.finalConditions
    delList = np.array([pc.pcIS.iG6P_c, pc.pcIS.iGLC_c, pc.pcIS.iFUM_c,
                        pc.pcIS.iFUM_i, pc.pcIS.iPYR_c, pc.pcIS.iK_c, pc.pcIS.iMg_c,
                        pc.pcIS.iH_c, pc.pcIS.iK_i, pc.pcIS.iMg_i, pc.pcIS.iH_i,
                        pc.pcIS.iCO2tot_x, pc.pcIS.iO2_x, pc.pcIS.iAMP_x], dtype=int)
    if not all(delList[i] >= delList[i+1] for i in range(len(delList)-1)):
        raise outOfOrder() ## Performs check for what follows that uses delList
    k = -1
    for i in range(len(y)): ## So that we do not have to provide every state
        ## Depends heavily on the indices being in decreasing order
        k = chooseK(i, k, delList)
        if y[i] != x[k]: ## If the state provided differs from the default final condition,
            ## prefer the value provided
            x[k] = y[i]
    result = equations.conservationEqs(x, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType)
    result = np.delete(result, delList) ## For the linearization, we only wish to provide
    ## the states that we've not excluded
    return result

def main():
    jac = nd.Jacobian(f)
    equi = np.array(pc.finalConditions)
    ## List of state variables we delete from final condition calculation
    delList = np.array([pc.pcIS.iG6P_c, pc.pcIS.iGLC_c, pc.pcIS.iFUM_c,
                        pc.pcIS.iFUM_i, pc.pcIS.iPYR_c, pc.pcIS.iK_c, pc.pcIS.iMg_c,
                        pc.pcIS.iH_c, pc.pcIS.iK_i, pc.pcIS.iMg_i, pc.pcIS.iH_i,
                        pc.pcIS.iCO2tot_x, pc.pcIS.iO2_x, pc.pcIS.iAMP_x], dtype = int)
    equi = np.delete(equi, delList)
    linearization = jac(equi)
    pseudo = ps.NonnormalMeshgrid(linearization,
                                  real_min=-5000000, real_max=5000000, real_n=1001,
                                  imag_min=-5000000, imag_max=5000000, imag_n=1001)
    pyplot.imshow(pseudo.Vals, extent=(-5000000, 5000000, -5000000, 5000000))
    pyplot.colorbar()
    pyplot.xlabel("Real Part")
    pyplot.ylabel("Imaginary Part")
    pyplot.show()

main()