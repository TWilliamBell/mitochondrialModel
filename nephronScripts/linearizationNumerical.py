import numpy as np
import numdifftools as nd
import sympy as sym

import equations
import pc

np.set_printoptions(threshold=10000)

J_AtC = 1.256e-3 ## Used by Edwards et al.
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

def f(y): ## Differential equations, with optional arguments specified
    return equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType)

def main():
    jac = nd.Jacobian(f)
    linearization = jac(pc.finalConditions)
    vals, vecs = np.linalg.eig(linearization)
    print(repr(linearization))
    #print(vals) ## Print the eigenvalues
    #print(vecs) ## Print the eigenvectors
    #charPol = np.poly(linearization) ## Routh-Hurwitz corroborates the
    ## existence of positive eigenvalues.
    #linearization = sym.Matrix(linearization)
    #null = np.array(linearization.nullspace())
    #for i in range(len(null)):
    #    print(np.matmul(np.array(linearization), null[i]))

main()