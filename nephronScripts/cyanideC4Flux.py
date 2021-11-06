import pandas as pd
import numpy as np

import fluxes as fl
import fluxesmTAL as flm
import pc as pc

def main():
    cyanide = pd.read_feather("../results/resultsCyanide.feather")
    baseline = pd.read_feather("../results/resultsATP.feather")
    baselinemTAL = pd.read_feather("../results/resultsmTAL.feather")
    cyanidemTAL = pd.read_feather("../results/resultsCyanidemTAL.feather")

    cyanide = np.array(cyanide.tail(1))[0]
    baseline = np.array(baseline.tail(1))[0]
    cyanidemTAL = np.array(cyanidemTAL.tail(1))[0]
    baselinemTAL = np.array(baselinemTAL.tail(1))[0]

    cyanide = np.delete(cyanide, 0)
    baseline = np.delete(baseline, 0)
    cyanidemTAL = np.delete(cyanidemTAL, 0)
    baselinemTAL = np.delete(baselinemTAL, 0)
    indHLe = [pc.pcIS.iH_i, pc.pcIS.iH_x, pc.pcIS.idPsi]
    # print(cyanide[0])
    # print(cyanide[indHLe])
    cyanideFlux = fl.fluxes(cyanide, pc.params, 1,
                                 w = [1., 1., 0.002, 1.])
    baselineFlux = fl.fluxes(baseline, pc.params, 1,
                                 w = [1., 1., 1., 1.])
    pc.params[38] = pc.params[38] / 7.0
    pc.pcPC.k_O2 = pc.pcPC.k_O2 / 10.0
    cyanidemTALFlux = flm.fluxesmTAL(cyanidemTAL, pc.params, 1,
                                 w = [1., 1., 0.002, 1.0])
    baselinemTALFlux = flm.fluxesmTAL(baselinemTAL, pc.params, 1,
                                 w = [1., 1., 1., 1.0])
    ## We investigate the Proton Potential Gradient-related fluxes.
    indPsi = np.array([pc.pcFl.J_C1, pc.pcFl.J_C3, pc.pcFl.J_C4, pc.pcFl.J_F1,
              pc.pcFl.J_ANT, pc.pcFl.J_Hle, pc.pcFl.J_ASP_GLU, pc.pcFl.J_Kle])

    print(baselineFlux[indPsi])
    print(cyanideFlux[indPsi])

    print(baselinemTALFlux[indPsi])
    print(cyanidemTALFlux[indPsi])

    indNADH = np.array([pc.pcFl.J_pdh, pc.pcFl.J_isod, pc.pcFl.J_akgd,
                        pc.pcFl.J_mdh])

    print(baselineFlux[indNADH])
    print(cyanideFlux[indNADH])

    print(baselinemTALFlux[indNADH])
    print(cyanidemTALFlux[indNADH])

    print("PO Ratio PT:")
    print(baselineFlux[pc.pcFl.J_F1]/baselineFlux[pc.pcFl.J_C4])
    print("PO Ratio mTAL:")
    print(baselinemTALFlux[pc.pcFl.J_F1]/baselinemTALFlux[pc.pcFl.J_C4])
    # print(cyanideFlux)#[pc.pcFl.J_Hle])
    # print(cyanideHLeakFlux)#[pc.pcFl.J_Hle])
    # print(cyanideFlux[indPsi]/baselineFlux[indPsi])
    ## We check the ratio of the fluxes between Cyanide poisoned case
    ## and the baseline, there are multiple that increase by orders of
    ## magnitude, but none that have large absolute value.
    # np.seterr(divide='ignore', invalid='ignore')
    # ratio = np.divide(cyanideFlux, baselineFlux)
    # print(ratio)
    # ratio = np.divide(cyanidemTALFlux, baselinemTALFlux)
    # print(ratio)

    # boole = ratio > 10
    # print(ratio)
    # print(boole)
    # ind = np.linspace(0, 45, 46)
    # print(ind[boole])
    # print(cyanideFlux[boole])
    # print(baselineFlux[boole])

main()