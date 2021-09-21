import numdifftools as nd
import numpy as np
from time import time

import equations
import pc

np.random.seed(1)

J_AtC = 1.256e-3 ## Used by Edwards et al.
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

def f(y): ## Differential equations, with optional arguments specified
    return equations.conservationEqs(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType)

def main():
    nIt = 1000
    eps = 1e-5
    jac = nd.Jacobian(f)
    results = np.zeros(shape = (nIt, pc.pc.N_states), dtype = np.complex)
    for i in range(nIt):
        while True: ## Generate strictly positive perturbations
            perturb = np.random.normal(0., 1.0, pc.pc.N_states) ## Guarantees uniformly
            ## distributed points over the hypersphere
            perturb = (1e-3/np.linalg.norm(perturb))*(perturb)
            for j in range(pc.pc.N_states):
                if pc.finalConditions[j] < eps:
                    perturb[j] = np.abs(perturb[j])
            point = pc.finalConditions+perturb
            if all(point[i] >= 0 for i in range(len(point))):
                break
        vals = np.linalg.eigvals(jac(point)) ## Take eigenvalues at perturbed points
        results[i, ] = vals
        print(i)
    np.savetxt("../results/perturbedEqui.csv", results, delimiter = ",")

main()
