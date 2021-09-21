import scipy.optimize as sci
import numpy as np
#import pyswarm as ps

import equations
import pc

#np.random.seed(1)

J_AtC = 1.256e-3 ## Used by Edwards et al.
ExpType = 1 ## in vivo = Pyruvate in cytoplasm clamped, cytoplasm has specified water
## volume
StateType = 1 ## Default, remaining Pyruvate concentrations not clamped

def f(y): ## Differential equations, with optional arguments specified
    return equations.conservationEqs1(y, J_AtC = J_AtC,
                              ExpType = ExpType,
                              StateType = StateType)

def main():
    ## Naive

    a = sci.root(f, x0 = pc.finalConditions, method = 'lm') ## lm, Krylov
    print(a.x)
    print(f(a.x))
    print(np.linalg.norm(f(a.x)))
    print(a.success)
    ## PSO
    #ub = pc.finalConditions
    #lb = np.zeros(63)
    #for i in range(len(ub)):
    #    ub[i] = np.max([ub[i]*1.1, 1e-3])
    # steadyState = ps.pso(cost,
    #                      lb = lb+1e-13,
    #                      ub = ub,
    #                      maxiter = 200,
    #                      swarmsize = 1000,
    #                      minstep = 0,
    #                      phig = 1.25,
    #                      omega = 1)

    #print(steadyState[1])
    #print(f(steadyState[0]))
    #print(np.linalg.norm(f(steadyState.x)/steadyState.x))

main()