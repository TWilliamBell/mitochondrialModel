import sklearn.cluster as sk
import numpy as np
import numdifftools as nd
import collections

import pc
import equations

J_AtC = 1.256e-3 ## Used by Edwards et al.
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped,
## cytoplasm has specified water volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

def f(y): ## Differential equations, with optional arguments specified
    return equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType)

## Finds the eigenvalues of the linearization at the steady state
## and cluster them using the DBSCAN algorithm.  The -1 values are
## unclustered, the 0 eigenvalues correspond to a repeated eigenvalue
## of roughly -3.98e03, the 1 eigenvalues are clustered around
## -2.5e-02, and the 2 eigenvalues are clustered around zero.
def main():
    jac = nd.Jacobian(f)
    linearization = jac(pc.finalConditions)
    vals, vecs = np.linalg.eig(linearization)
    clust = np.array([np.real(vals), np.imag(vals)]).transpose()
    clustered = sk.DBSCAN(eps = 40, min_samples = 2).fit(clust)
    print(vals)
    print(clustered.labels_)
    table = collections.Counter(clustered.labels_)
    print(table)

main()
