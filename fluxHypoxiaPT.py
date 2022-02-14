import scipy.integrate as sci
import numpy as np
import pandas as pd

import equations
import fluxes as fl
import pc

J_AtC = 1.256e-3

class notFound(Exception):
    def __init__(self):
        self.message = "Not found"

def ret_index(t, crit = 12000):
    for i in range(len(t)):
        if t[i] > crit:
            return i
    raise notFound


def main():
    soln = pd.read_csv("../results/resultsHypoxia1.csv")

    for i in range(5, 35):
        crit = i*1000+1000
        crucial_index = ret_index(soln['t'].to_numpy(), crit = crit)
        state_time = np.delete(soln.loc[crucial_index, ].to_numpy(), [0, 1])
        flux_time = fl.fluxes(state_time, param = pc.params, ExpType = 1)
        #print(flux_time)

        Rm_cyto = 0.25/0.72

        ATPfl = [-Rm_cyto*flux_time[pc.pcFl.J_ATP],
             flux_time[pc.pcFl.J_CKe],
             flux_time[pc.pcFl.J_AKe],
             -(1.63*J_AtC)*(state_time[pc.pcIS.iATP_c]/(1.486e-3+state_time[pc.pcIS.iATP_c]))]
        ## See what the balance is between the main driver of ATP flux
        print(ATPfl[0]+ATPfl[3])

main()