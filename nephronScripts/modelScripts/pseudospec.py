import numdifftools as nd
import pseudopy as ps
from matplotlib import pyplot
from scipy.linalg import eigvals
import pandas as pd

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

def main():
    jac = nd.Jacobian(f)
    linearization = jac(pc.finalConditions)

    ## Pseudospectrum
    pseudo = ps.NonnormalMeshgrid(linearization,
                                  real_min = -25, real_max = 25, real_n = 101,
                                  imag_min = -20, imag_max = 20, imag_n = 101)
    pyplot.imshow(pseudo.Vals, extent = (-25, 25, -20, 20))
    pyplot.colorbar()
    pyplot.xlabel("Real Part")
    pyplot.ylabel("Imaginary Part")
    pyplot.show()

    ## Save linearization
    linearization = pd.DataFrame(linearization)
    linearization.to_csv("../results/linearization.csv")



main()

