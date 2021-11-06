import numpy as np
import numdifftools as nd
import sympy as sym

import equations
import pc

J_AtC = 1.256e-3 ## Used by Edwards et al.
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

def f(y): ## Differential equations, with optional arguments specified
    return equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType)


def dL():
    return np.array([pc.pcIS.iG6P_c, ## Confirmed directly
                     pc.pcIS.iGLC_c, ## Confirmed directly
                     pc.pcIS.iFUM_c, ## Confirmed directly
                     pc.pcIS.iFUM_i, ## Confirmed directly
                     pc.pcIS.iK_i, ## Confirmed directly
                     pc.pcIS.iMg_i, ## Confirmed directly
                     pc.pcIS.iH_i ## Confirmed directly
                     ], dtype = int)

def main():
    jac = nd.Jacobian(f)
    linearization = jac(pc.finalConditions)
    vals, vecs = np.linalg.eig(linearization)
    for i in range(62):
        eigvec = np.zeros((63, 1))
        eigvec[i] = 1
        eigvec[i+1] = 1
        Ax = np.matmul(linearization, eigvec)
        if np.linalg.norm(Ax) < 1e-10:
            print(i)
            print(Ax)
main()