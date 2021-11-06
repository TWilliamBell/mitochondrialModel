import numpy as np
import time as time

from classDefs import negativeConcentration, timeError
from pc import *
import fluxes as fl
import fluxesmTAL as flt

def conservationEqs(x, J_AtC, ExpType = 1, StateType = 1,
                    w = [1., 1., 1., 1.], timeStart = 1,
                    capacitanceWeight = 1., potassiumW = 1.,
                    DCA = 1., tubule = "PT", glyc = 0.0,
                    cana = 1., params = params):
    """The conservation equations first computes the fluxes, and then
    calculates the differential equations based on the fluxes returned."""

    if (not timeStart == 1) and np.abs(timeStart-time.time()) > 3600:
        raise timeError()

    for i in range(len(x)):
        if x[i] < 0:
            x[i] = 1e-20
    if tubule == "PT":
        Jsol1 = fl.fluxes(x, params, ExpType = ExpType, w = w,
                      potassiumW = potassiumW, DCA = DCA, cana = 1.)
    elif tubule == "mTAL":
        pcPC.W_x = 0.58*0.9
        pcPC.W_i = 0.58*0.1
        pcPC.W_m = 0.58
        Jsol1 = flt.fluxesmTAL(x, params, ExpType = ExpType, w = w,
                      potassiumW = potassiumW, DCA = DCA)

    ## Experimentally-dependent parameters
    if tubule == "PT":
        if ExpType == 1:
            W_c = 0.77 # cytoplasm water fraction - LeFurgey et al. 1991
            Vmito = 0.250 # (ml mito / ml cell)
            Vcyto = 0.720 # (ml cyto / ml cell)
            Rm_cyto = Vmito / Vcyto # Volume ratio mito volume / cytoplasm volume
            Rm_cell = Vmito # Volume ratio mito volume / cell volume
            Rc_cell = Vcyto # Volume ratio of cytosol / cell volume
        elif ExpType == 3:
            W_c = 0.77 # cytoplasm water fraction - LeFurgey et al. 1991
            Vmito = 0.313 # (ml mito / ml cell)
            Vcyto = 0.658 # (ml cyto / ml cell)
            Rm_cyto = Vmito / Vcyto # Volume ratio mito volume / cytoplasm volume
            Rm_cell = Vmito # Volume ratio mito volume / cell volume
            Rc_cell = Vcyto # Volume ratio of cytosol / cell volume
        elif ExpType == 2: # Bose et al.'s experiments (JBC, 2003)
            W_c = 80 * 10 # Buffer water space
        else:
            W_c = 80
    elif tubule == "mTAL":
        if ExpType == 1:
            W_c = 0.77 # cytoplasm water fraction - LeFurgey et al. 1991
            Vmito = 0.32 # (ml mito / ml cell)
            Vcyto = 0.59 # (ml cyto / ml cell)
            Rm_cyto = Vmito / Vcyto # Volume ratio mito volume / cytoplasm volume
            Rm_cell = Vmito # Volume ratio mito volume / cell volume
            Rc_cell = Vcyto # Volume ratio of cytosol / cell volume
        elif ExpType == 2: # Bose et al.'s experiments (JBC, 2003)
            W_c = 80 * 10 # Buffer water space
        else:
            W_c = 80

    ## Computing time derivatives of state variables

    f = np.zeros(pc.N_states)

    ## (i) Matrix species and dPsi
    f[pcIS.idPsi] = 1.48e5*(4*Jsol1[pcFl.J_C1] + 2*Jsol1[pcFl.J_C3] +
        4*Jsol1[pcFl.J_C4] - pcPC.n_A*Jsol1[pcFl.J_F1] - Jsol1[pcFl.J_ANT]
        - Jsol1[pcFl.J_Hle] - Jsol1[pcFl.J_Kle] + Jsol1[pcFl.J_ASP_GLU])/\
        capacitanceWeight
    f[pcIS.iATP_x] = (Jsol1[pcFl.J_ndk] + Jsol1[pcFl.J_F1] - Jsol1[pcFl.J_ANT])/pcPC.W_x
    f[pcIS.iADP_x]  = (-Jsol1[pcFl.J_ndk] - Jsol1[pcFl.J_F1] + Jsol1[pcFl.J_ANT])/pcPC.W_x
    f[pcIS.iAMP_x]  = 0./pcPC.W_x
    f[pcIS.iGTP_x]  = (+Jsol1[pcFl.J_scoas] - Jsol1[pcFl.J_ndk])/pcPC.W_x
    f[pcIS.iGDP_x]  = (-Jsol1[pcFl.J_scoas] + Jsol1[pcFl.J_ndk])/pcPC.W_x
    f[pcIS.iPi_x]   = (-Jsol1[pcFl.J_scoas] - Jsol1[pcFl.J_F1] + Jsol1[pcFl.J_Pi1]
        - Jsol1[pcFl.J_MAL_PI])/pcPC.W_x
    f[pcIS.iNADH_x] = (+Jsol1[pcFl.J_pdh] + Jsol1[pcFl.J_isod] + Jsol1[pcFl.J_akgd] +
        Jsol1[pcFl.J_mdh] - Jsol1[pcFl.J_C1])/pcPC.W_x
    f[pcIS.iQH2_x]  = (+Jsol1[pcFl.J_sdh] + Jsol1[pcFl.J_C1] - Jsol1[pcFl.J_C3])/pcPC.W_x
    f[pcIS.iPYR_x]  = (-Jsol1[pcFl.J_pdh] + Jsol1[pcFl.J_PYR_H])/pcPC.W_x
    f[pcIS.iACCOA_x] = (-Jsol1[pcFl.J_cits] + Jsol1[pcFl.J_pdh])/pcPC.W_x
    f[pcIS.iCIT_x]  = (+Jsol1[pcFl.J_cits] - Jsol1[pcFl.J_acon] + Jsol1[pcFl.J_CIT_MAL])/ \
        pcPC.W_x
    f[pcIS.iICIT_x] = (+Jsol1[pcFl.J_acon] - Jsol1[pcFl.J_isod])/pcPC.W_x
    f[pcIS.iAKG_x]  = (+Jsol1[pcFl.J_isod] - Jsol1[pcFl.J_akgd] - Jsol1[pcFl.J_got] +
        Jsol1[pcFl.J_AKG_MAL])/pcPC.W_x
    f[pcIS.iSCOA_x] = (+Jsol1[pcFl.J_akgd] - Jsol1[pcFl.J_scoas])/pcPC.W_x
    f[pcIS.iCOASH_x]  = (-Jsol1[pcFl.J_pdh] - Jsol1[pcFl.J_akgd] + Jsol1[pcFl.J_scoas]
        + Jsol1[pcFl.J_cits])/pcPC.W_x
    f[pcIS.iSUC_x]  = (+Jsol1[pcFl.J_scoas] - Jsol1[pcFl.J_sdh] +
        Jsol1[pcFl.J_SUC_MAL])/pcPC.W_x
    f[pcIS.iFUM_x]  = (+Jsol1[pcFl.J_sdh] - Jsol1[pcFl.J_fum])/pcPC.W_x
    f[pcIS.iMAL_x]  = (+Jsol1[pcFl.J_fum] - Jsol1[pcFl.J_mdh] + Jsol1[pcFl.J_MAL_PI] -
        Jsol1[pcFl.J_AKG_MAL] - Jsol1[pcFl.J_CIT_MAL] - Jsol1[pcFl.J_SUC_MAL])/pcPC.W_x
    f[pcIS.iOAA_x]  = (-Jsol1[pcFl.J_cits] + Jsol1[pcFl.J_mdh] + Jsol1[pcFl.J_got])/pcPC.W_x

    f[pcIS.iGLU_x]  = (+Jsol1[pcFl.J_got] + Jsol1[pcFl.J_GLU_H]
        - Jsol1[pcFl.J_ASP_GLU])/pcPC.W_x
    f[pcIS.iASP_x]  = (-Jsol1[pcFl.J_got] + Jsol1[pcFl.J_ASP_GLU])/pcPC.W_x

    ## (ii) IM space species
    f[pcIS.iCred_i] = (+2*Jsol1[pcFl.J_C3] - 2*Jsol1[pcFl.J_C4])/pcPC.W_i
    f[pcIS.iATP_i] = (Jsol1[pcFl.J_ATP] + Jsol1[pcFl.J_ANT] + Jsol1[pcFl.J_AKi])/pcPC.W_i
    f[pcIS.iADP_i] = (Jsol1[pcFl.J_ADP] - Jsol1[pcFl.J_ANT] - 2*Jsol1[pcFl.J_AKi])/pcPC.W_i
    f[pcIS.iAMP_i] = (Jsol1[pcFl.J_AMP] + Jsol1[pcFl.J_AKi])/pcPC.W_i
    f[pcIS.iPi_i]  = (-Jsol1[pcFl.J_Pi1] + Jsol1[pcFl.J_Pi2] +
        Jsol1[pcFl.J_MAL_PI])/pcPC.W_i
    f[pcIS.iPYR_i] = (-Jsol1[pcFl.J_PYR_H] + Jsol1[pcFl.J_PYRt])/pcPC.W_i
    f[pcIS.iCIT_i] = (-Jsol1[pcFl.J_CIT_MAL] + Jsol1[pcFl.J_CITt])/pcPC.W_i
    f[pcIS.iICIT_i] = (Jsol1[pcFl.J_ICITt])/pcPC.W_i
    f[pcIS.iAKG_i] = (-Jsol1[pcFl.J_AKG_MAL] + Jsol1[pcFl.J_AKGt])/pcPC.W_i
    f[pcIS.iSUC_i] = (+ Jsol1[pcFl.J_SUCt] - Jsol1[pcFl.J_SUC_MAL])/pcPC.W_i
    f[pcIS.iFUM_i] = (Jsol1[pcFl.J_FUMt])/pcPC.W_i
    f[pcIS.iMAL_i] = (-Jsol1[pcFl.J_MAL_PI] + Jsol1[pcFl.J_MALt] +
        Jsol1[pcFl.J_AKG_MAL] + Jsol1[pcFl.J_CIT_MAL] +
        Jsol1[pcFl.J_SUC_MAL])/pcPC.W_i
    f[pcIS.iGLU_i] = (-Jsol1[pcFl.J_GLU_H] + Jsol1[pcFl.J_ASP_GLU] +
        Jsol1[pcFl.J_GLUt])/pcPC.W_i
    f[pcIS.iASP_i] = (-Jsol1[pcFl.J_ASP_GLU] + Jsol1[pcFl.J_ASPt])/pcPC.W_i

    ## (iii) Buffer species
    #  d[ADP]/dt, d[ATP]/dt and d[PI]/dt in cytoplasm are assumed to be
    #  consumed by Hexokinase reaction
    #  W_c   [=] buffer water/mito || cyto water/cyto
    #  J_ATP [=] M/s/l mito
    #  J_HK  [=] M/s/l buffer
    #  J_AtC, J_CKe [=] M/s/l cyto
    #  f [=] M/s/l cyto water
    if (ExpType == 2) or (ExpType == 5) or (ExpType == 6):
        f[pcIS.iPYR_c] = (-Jsol1[pcFl.J_PYRt])/W_c
        f[pcIS.iCIT_c] = (-Jsol1[pcFl.J_CITt])/W_c
        f[pcIS.iICIT_c] = (-Jsol1[pcFl.J_ICITt])/W_c
        f[pcIS.iAKG_c] = (-Jsol1[pcFl.J_AKGt])/W_c
        f[pcIS.iSUC_c] = (-Jsol1[pcFl.J_SUCt])/W_c
        f[pcIS.iFUM_c] = (-Jsol1[pcFl.J_FUMt])/W_c
        f[pcIS.iMAL_c] = (-Jsol1[pcFl.J_MALt])/W_c
        f[pcIS.iGLU_c] = (-Jsol1[pcFl.J_GLUt])/W_c
        f[pcIS.iASP_c] = (-Jsol1[pcFl.J_ASPt])/W_c
        f[pcIS.iPi_c]  = (-Jsol1[pcFl.J_Pi2])/W_c
        f[pcIS.iATP_c] = -Jsol1[pcFl.J_HK]/(W_c/(1+W_c))- \
            Jsol1[pcFl.J_ATP]/W_c # SET TO ZERO IN P/O MEASUREMENTS
        f[pcIS.iADP_c] = +Jsol1[pcFl.J_HK]/(W_c/(1+W_c))- \
            Jsol1[pcFl.J_ADP]/W_c # SET TO ZERO IN P/O MEASUREMENTS
        f[pcIS.iGLC_c] = -Jsol1[pcFl.J_HK]/(W_c/(1+W_c))
        f[pcIS.iG6P_c] = +Jsol1[pcFl.J_HK]/(W_c/(1+W_c))
    else:
        f[pcIS.iPYR_c] = (-Rm_cyto*Jsol1[pcFl.J_PYRt])/W_c
        f[pcIS.iCIT_c] = (-Rm_cyto*Jsol1[pcFl.J_CITt])/W_c
        f[pcIS.iICIT_c] = (-Rm_cyto*Jsol1[pcFl.J_ICITt])/W_c
        f[pcIS.iAKG_c] = (-Rm_cyto*Jsol1[pcFl.J_AKGt])/W_c
        f[pcIS.iSUC_c] = (-Rm_cyto*Jsol1[pcFl.J_SUCt])/W_c
        f[pcIS.iFUM_c] = (-Rm_cyto*Jsol1[pcFl.J_FUMt])/W_c
        f[pcIS.iMAL_c] = (-Rm_cyto*Jsol1[pcFl.J_MALt])/W_c
        f[pcIS.iGLU_c] = (-Rm_cyto*Jsol1[pcFl.J_GLUt])/W_c
        f[pcIS.iASP_c] = (-Rm_cyto*Jsol1[pcFl.J_ASPt])/W_c
        f[pcIS.iPi_c]  = (-Rm_cyto*Jsol1[pcFl.J_Pi2]+J_AtC-glyc)/W_c
        f[pcIS.iATP_c] = (-Rm_cyto*Jsol1[pcFl.J_ATP]-J_AtC+glyc+
            Jsol1[pcFl.J_CKe]+Jsol1[pcFl.J_AKe])/W_c
        f[pcIS.iADP_c] = (-Rm_cyto*Jsol1[pcFl.J_ADP]+J_AtC-glyc-
            Jsol1[pcFl.J_CKe]-2*Jsol1[pcFl.J_AKe])/W_c
        f[pcIS.iPCr_c] = (-Jsol1[pcFl.J_CKe])/W_c
        f[pcIS.iAMP_c] = (-Rm_cyto*Jsol1[pcFl.J_AMP]+Jsol1[pcFl.J_AKe])/W_c

    if (StateType == 0):
        f[pcIS.iPYR_c] = 0.
        f[pcIS.iMAL_c] = 0.
        f[pcIS.iCIT_c] = 0.
        f[pcIS.iICIT_c] = 0.
        f[pcIS.iAKG_c] = 0.
        f[pcIS.iSUC_c] = 0.
        f[pcIS.iFUM_c] = 0.

    if (ExpType == 1) or (ExpType == 3) or (ExpType == 4):
        f[pcIS.iPYR_c] = 0.
    elif (ExpType == 5): # for uncoupling
        f[pcIS.iPYR_c] = 0.
        f[pcIS.iMAL_c] = 0.
        f[pcIS.iCIT_c] = 0.
        f[pcIS.iICIT_c] = 0.
        f[pcIS.iAKG_c] = 0.
        f[pcIS.iSUC_c] = 0.
        f[pcIS.iFUM_c] = 0.

    # %% Computing dH/dt, dMg/dt, and dK/dt
    # % Calculate dH/dt, dMg/dt, dK/dt according to Dan's strategy (available
    # % in the chapter "Biochemical Reaction Networks" in Beard and Qian's book)
    #
    # % assume [H+], [Mg2+], [K+] constant in IM and cytoplasm/buffer space
    # % all the below calculation if for dH_x/dt, dMg_x/dt, dK_x/dt
    #
    # % take concentration of reactants located in the mito. matrix
    # % from state variables

    xx = np.zeros(pc.N_reactant)
    xx[pcIR.iH] = x[pcIS.iH_x]
    xx[pcIR.iATP] = x[pcIS.iATP_x]
    xx[pcIR.iADP] = x[pcIS.iADP_x]
    xx[pcIR.iAMP] = x[pcIS.iAMP_x]
    xx[pcIR.iGTP] = x[pcIS.iGTP_x]
    xx[pcIR.iGDP] = x[pcIS.iGDP_x]
    xx[pcIR.iPi] = x[pcIS.iPi_x]
    xx[pcIR.iNADH] = x[pcIS.iNADH_x]
    if tubule == "mTAL":
        xx[pcIR.iNAD] = pcPC.NADtotmTAL - x[pcIS.iNADH_x]
    else:
        xx[pcIR.iNAD] = pcPC.NADtot - x[pcIS.iNADH_x]
    xx[pcIR.iQH2] = x[pcIS.iQH2_x]
    xx[pcIR.iCOQ] = pcPC.Qtot - x[pcIS.iQH2_x]
    xx[pcIR.iPYR] = x[pcIS.iPYR_x]
    xx[pcIR.iOAA] = x[pcIS.iOAA_x]
    xx[pcIR.iACCOA] = x[pcIS.iACCOA_x]
    xx[pcIR.iCIT] = x[pcIS.iCIT_x]
    xx[pcIR.iICIT] = x[pcIS.iICIT_x]
    xx[pcIR.iAKG] = x[pcIS.iAKG_x]
    xx[pcIR.iSCOA] = x[pcIS.iSCOA_x]
    xx[pcIR.iCOASH] = x[pcIS.iCOASH_x]
    xx[pcIR.iSUC] = x[pcIS.iSUC_x]
    xx[pcIR.iFUM] = x[pcIS.iFUM_x]
    xx[pcIR.iMAL] = x[pcIS.iMAL_x]
    xx[pcIR.iGLU] = x[pcIS.iGLU_x]
    xx[pcIR.iASP] = x[pcIS.iASP_x]
    xx[pcIR.iK] = x[pcIS.iK_x]
    xx[pcIR.iMg] = x[pcIS.iMg_x]
    xx[pcIR.iO2] = x[pcIS.iO2_x]
    xx[pcIR.iFADH2] = pcPC.FADtot/2
    xx[pcIR.iFAD] = pcPC.FADtot - xx[pcIR.iFADH2]
    xx[pcIR.iCO2tot] = x[pcIS.iCO2tot_x]

    dxxdt = np.zeros(pc.N_reactant)
    # dxxdt[pcIR.iH] = f[pcIS.iH_x]
    dxxdt[pcIR.iATP] = f[pcIS.iATP_x]
    dxxdt[pcIR.iADP] = f[pcIS.iADP_x]
    dxxdt[pcIR.iAMP] = f[pcIS.iAMP_x]
    dxxdt[pcIR.iGTP] = f[pcIS.iGTP_x]
    dxxdt[pcIR.iGDP] = f[pcIS.iGDP_x]
    dxxdt[pcIR.iPi] = f[pcIS.iPi_x]
    dxxdt[pcIR.iNADH] = f[pcIS.iNADH_x]
    dxxdt[pcIR.iNAD] = -f[pcIS.iNADH_x]
    dxxdt[pcIR.iQH2] = f[pcIS.iQH2_x]
    dxxdt[pcIR.iCOQ] = - f[pcIS.iQH2_x]
    dxxdt[pcIR.iPYR] = f[pcIS.iPYR_x]
    dxxdt[pcIR.iOAA] = f[pcIS.iOAA_x]
    dxxdt[pcIR.iACCOA] = f[pcIS.iACCOA_x]
    dxxdt[pcIR.iCIT] = f[pcIS.iCIT_x]
    dxxdt[pcIR.iICIT] = f[pcIS.iICIT_x]
    dxxdt[pcIR.iAKG] = f[pcIS.iAKG_x]
    dxxdt[pcIR.iSCOA] = f[pcIS.iSCOA_x]
    dxxdt[pcIR.iCOASH] = f[pcIS.iCOASH_x]
    dxxdt[pcIR.iSUC] = f[pcIS.iSUC_x]
    dxxdt[pcIR.iFUM] = f[pcIS.iFUM_x]
    dxxdt[pcIR.iMAL] = f[pcIS.iMAL_x]
    dxxdt[pcIR.iGLU] = f[pcIS.iGLU_x]
    dxxdt[pcIR.iASP] = f[pcIS.iASP_x]
    # dxxdt[pcIR.iK]   = f[pcIS.iK_x]
    # dxxdt[pcIR.iMg]  = f[pcIS.iMg_x]
    dxxdt[pcIR.iCred] = f[pcIS.iCred_i]
    dxxdt[pcIR.iCox] = -f[pcIS.iCred_i]
    dxxdt[pcIR.iO2]  = f[pcIS.iO2_x]
    # dxxdt[pcIR.iH2O] = 0
    # dxxdt[pcIR.iFADH2] = 0
    dxxdt[pcIR.iFAD] = -dxxdt[pcIR.iFADH2]
    dxxdt[pcIR.iCO2tot] = f[pcIS.iCO2tot_x]

    # Necessary concentrations
    H_x = x[pcIS.iH_x]
    K_x = x[pcIS.iK_x]
    Mg_x = x[pcIS.iMg_x]

    ## K values for reference species
    ## Repeated code from fluxes.py
    Kh = [float('inf')] * pc.N_reactant
    Km = [float('inf')] * pc.N_reactant
    Kk = [float('inf')] * pc.N_reactant
    Kh[pcIR.iATP] = 10 ** (-6.59)
    Km[pcIR.iATP] = 10 ** (-3.82)
    Kk[pcIR.iATP] = 10 ** (-1.013)
    Kh[pcIR.iADP] = 10 ** (-6.42)
    Km[pcIR.iADP] = 10 ** (-2.79)
    Kk[pcIR.iADP] = 10 ** (-0.882)
    Kh[pcIR.iAMP] = 10 ** (-6.22)
    Km[pcIR.iAMP] = 10 ** (-1.86)
    Kk[pcIR.iAMP] = 10 ** (-0.6215)
    Kh[pcIR.iGTP] = Kh[pcIR.iATP]
    Km[pcIR.iGTP] = Km[pcIR.iATP]
    Kk[pcIR.iGTP] = Kk[pcIR.iATP]
    Kh[pcIR.iGDP] = Kh[pcIR.iADP]
    Km[pcIR.iGDP] = Km[pcIR.iADP]
    Kk[pcIR.iGDP] = Kk[pcIR.iADP]
    Kh[pcIR.iPi] = 10 ** (-6.71)
    Km[pcIR.iPi] = 10 ** (-1.69)
    Kk[pcIR.iPi] = 10 ** (+0.0074)
    Kh[pcIR.iCOASH] = 10 ** (-8.13)
    Km[pcIR.iOAA] = 10 ** (-0.8629)
    Kh[pcIR.iCIT] = 10 ** (-5.63)
    Km[pcIR.iCIT] = 10 ** (-3.37)
    Kk[pcIR.iCIT] = 10 ** (-0.339)
    Kh[pcIR.iICIT] = 10 ** (-5.64)
    Km[pcIR.iICIT] = 10 ** (-2.46)
    Kh[pcIR.iSCOA] = 10 ** (-3.96)
    Kh[pcIR.iSUC] = 10 ** (-5.13)
    Km[pcIR.iSUC] = 10 ** (-1.17)
    Kk[pcIR.iSUC] = 10 ** (-0.3525)
    Kh[pcIR.iFUM] = 10 ** (-4.10)
    Kh[pcIR.iMAL] = 10 ** (-4.75)
    Km[pcIR.iMAL] = 10 ** (-1.55)
    Kk[pcIR.iMAL] = 10 ** (+0.107)
    Kh[pcIR.iCO2tot] = 10 ** (-9.82)
    Km[pcIR.iPYR] = 10 ** (-1.02)
    Kh[pcIR.iG6P] = 10 ** (-5.91)
    # Kh(iGLU) = 10**(-4.25) # from Nelson & Cox, "Lehninger's Princinples of Biochemistry"
    # , p78
    Kh[pcIR.iGLU] = 10 ** (-4.06)  # 37 C, I = 0.15
    Km[pcIR.iGLU] = 10 ** (-1.82)
    # Kh(iASP) = 10 ^ (-3.65) # from Nelson & Cox, "Lehninger's Princinples of Biochemistry"
    # , p78
    Kh[pcIR.iASP] = 10 ** (-3.65)  # 37 C, I = 0.15
    Km[pcIR.iASP] = 10 ** (-2.32)

    ## # compute binding polynomials for reactants, also repeated code from fluxes.py
    P_x = [1.]*pc.N_reactant
    P_x[pcIR.iATP] = 1 + H_x/Kh[pcIR.iATP] + Mg_x/Km[pcIR.iATP] + K_x/Kk[pcIR.iATP]
    P_x[pcIR.iADP] = 1 + H_x/Kh[pcIR.iADP] + Mg_x/Km[pcIR.iADP] + K_x/Kk[pcIR.iADP]
    P_x[pcIR.iAMP] = 1 + H_x/Kh[pcIR.iAMP] + Mg_x/Km[pcIR.iAMP] + K_x/Kk[pcIR.iAMP]
    P_x[pcIR.iGTP] = P_x[pcIR.iATP]
    P_x[pcIR.iGDP] = P_x[pcIR.iADP]
    P_x[pcIR.iPi] = 1 + H_x/Kh[pcIR.iPi] + Mg_x/Km[pcIR.iPi] + K_x/Kk[pcIR.iPi]
    P_x[pcIR.iCOASH] = 1 + H_x/Kh[pcIR.iCOASH]
    P_x[pcIR.iOAA] = 1 + Mg_x/Km[pcIR.iOAA]
    P_x[pcIR.iCIT] = 1 + H_x/Kh[pcIR.iCIT] + Mg_x/Km[pcIR.iCIT] + K_x/Kk[pcIR.iCIT]
    P_x[pcIR.iICIT] = 1 + H_x/Kh[pcIR.iICIT] + Mg_x/Km[pcIR.iICIT]
    P_x[pcIR.iSCOA] = 1 + H_x/Kh[pcIR.iSCOA]
    P_x[pcIR.iSUC] = 1 + H_x/Kh[pcIR.iSUC] + Mg_x/Km[pcIR.iSUC] + K_x/Kk[pcIR.iSUC]
    P_x[pcIR.iFUM] = 1 + H_x/Kh[pcIR.iFUM]
    P_x[pcIR.iMAL] = 1 + H_x/Kh[pcIR.iMAL] + Mg_x/Km[pcIR.iMAL] + K_x/Kk[pcIR.iMAL]
    P_x[pcIR.iCO2tot] = 1 + H_x/Kh[pcIR.iCO2tot]
    P_x[pcIR.iPYR] = 1 + Mg_x/Km[pcIR.iPYR]
    P_x[pcIR.iG6P] = 1 + H_x/Kh[pcIR.iG6P]
    P_x[pcIR.iGLU] = 1 + H_x/Kh[pcIR.iGLU] + Mg_x/Km[pcIR.iGLU]
    P_x[pcIR.iASP] = 1 + H_x/Kh[pcIR.iASP] + Mg_x/Km[pcIR.iASP]

    P_x = np.array(P_x)
    H_x = np.array(H_x)
    Mg_x = np.array(Mg_x)
    K_x = np.array(K_x)
    Kh = np.array(Kh)
    Km = np.array(Km)
    Kk = np.array(Kk)

    # Partial Derivatives
    pHBpM = - np.sum( (H_x*xx/Kh)/(Km*np.power(P_x,2)) )
    pHBpK = - np.sum( (H_x*xx/Kh)/(Kk*np.power(P_x,2)) )
    pHBpH = np.sum( (1+Mg_x/Km+K_x/Kk)*xx/(Kh*np.power(P_x,2)) )
    pMBpH = - np.sum( (Mg_x*xx/Km)/(Kh*np.power(P_x,2)) )
    pMBpK = - np.sum( (Mg_x*xx/Km)/(Kk*np.power(P_x,2)) )
    pMBpM = np.sum( (1+H_x/Kh+K_x/Kk)*xx/(Km*np.power(P_x,2)) )
    pKBpH = - np.sum( (K_x*xx/Kk)/(Kh*np.power(P_x,2)) )
    pKBpM = - np.sum( (K_x*xx/Kk)/(Km*np.power(P_x,2)) )
    pKBpK = np.sum( (1+H_x/Kh+Mg_x/Km)*xx/(Kk*np.power(P_x,2)) )

    # PHI's:
    Phi_H = - np.sum( H_x*dxxdt/(Kh*P_x) ) + 1*(-1*Jsol1[pcFl.J_pdh] +
        2*Jsol1[pcFl.J_cits] + (-1)*Jsol1[pcFl.J_akgd] + Jsol1[pcFl.J_scoas]
        + Jsol1[pcFl.J_mdh]+ 1*(Jsol1[pcFl.J_PYR_H] + Jsol1[pcFl.J_GLU_H] +
        Jsol1[pcFl.J_CIT_MAL] - Jsol1[pcFl.J_ASP_GLU])-(4+1)*Jsol1[pcFl.J_C1] -
        (4-2)*Jsol1[pcFl.J_C3] - (2+2)*Jsol1[pcFl.J_C4] + (pcPC.n_A-1)*Jsol1[pcFl.J_F1]
        + 2*Jsol1[pcFl.J_Pi1] + 1*Jsol1[pcFl.J_Hle] - 1*Jsol1[pcFl.J_KH])/pcPC.W_x
    Phi_M = -np.sum( Mg_x*dxxdt/(Km*P_x) )
    Phi_K = - np.sum( K_x*dxxdt/(Kk*P_x) ) + 1*Jsol1[pcFl.J_Kle]/pcPC.W_x \
        + 1*Jsol1[pcFl.J_KH]/pcPC.W_x

    # alpha's:
    aH = 1 + pHBpH
    aM = 1 + pMBpM
    aK = 1 + pKBpK

    # add additional buffer for [H+]
    BX = 0.02 # M
    K_BX = 1e-7 # M
    aH = 1 + pHBpH + BX/K_BX/(1+H_x/K_BX)**2

    #  Denominator:
    D = aH*pKBpM*pMBpK + aK*pHBpM*pMBpH + aM*pHBpK*pKBpH - \
        aM*aK*aH - pHBpK*pKBpM*pMBpH - pHBpM*pMBpK*pKBpH

    # Derivatives for H, Mg, K:
    dH_xdt = ((pKBpM * pMBpK - aM * aK) * Phi_H +
        (aK * pHBpM - pHBpK * pKBpM) * Phi_M +
        (aM * pHBpK - pHBpM * pMBpK) * Phi_K) / D

    dMg_xdt = ((aK * pMBpH - pKBpH * pMBpK) * Phi_H +
        (pKBpH * pHBpK - aH * aK) * Phi_M +
        (aH * pMBpK - pHBpK * pMBpH) * Phi_K) / D

    dK_xdt = ((aM * pKBpH - pKBpM * pMBpH) * Phi_H +
        (aH * pKBpM - pKBpH * pHBpM) * Phi_M +
        (pMBpH * pHBpM - aH * aM) * Phi_K) / D

    f[pcIS.iH_x] = dH_xdt
    f[pcIS.iMg_x] = dMg_xdt
    f[pcIS.iK_x] = dK_xdt

    return f

def conservationEqs1(x, J_AtC, ExpType = 1, StateType = 1,
                    w = [1., 1., 1., 1.], timeStart = 1,
                    capacitanceWeight = 1.,
                    potassiumW = 1., tubule = "PT",
                    glyc = 0.0, DCA = 1., cana = 1.,
                    params = params):
    """The conservation equations first computes the fluxes, and then
    calculates the differential equations based on the fluxes returned.
    The function here as opposed to the previous version of it
    provides an adaptive response of ATP consumption to
    low ATP levels."""

    if timeStart == 1:
        timeStart = time.time()

    if np.abs(timeStart-time.time()) > 3600:
        raise timeError()

    for i in range(len(x)):
        if x[i] < 0:
            x[i] = 1e-20

    if tubule == "PT":
        Jsol1 = fl.fluxes(x, params, ExpType = ExpType, w = w,
                      potassiumW = potassiumW, DCA = DCA, cana = cana)
    elif tubule == "mTAL":
        pcPC.W_x = 0.58*0.9
        pcPC.W_i = 0.58*0.1
        pcPC.W_m = 0.58
        Jsol1 = flt.fluxesmTAL(x, params, ExpType = ExpType, w = w,
                      potassiumW = potassiumW, DCA = DCA)

    J_AtC = (1.63*J_AtC)*(x[pcIS.iATP_c]/(1.486e-3+x[pcIS.iATP_c]))
    if tubule == "mTAL": ## Glycolysis term - still not sure about it
        J_AtC = J_AtC-glyc #1.83e-5
    ## Changed to maintain postive invariance of ATP_c, which appears to be a
    ## source of problems for the model.  This should not change the normal
    ## behaviour, params chosen to guarantee that, but it means that for small
    ## concentrations the rate of [ATP]_c use is reduced.  K_m is that for
    ## Na-K-ATPase and the V_max is chosen so that the J_AtC at the ordinary
    ## concentration range is similar to predictions from nephron ATP
    ## consumption models.

    ## Experimentally-dependent parameters
    if tubule == "PT":
        if ExpType == 1:
            W_c = 0.77 # cytoplasm water fraction - LeFurgey et al. 1991
            Vmito = 0.250 # (ml mito / ml cell)
            Vcyto = 0.720 # (ml cyto / ml cell)
            Rm_cyto = Vmito / Vcyto # Volume ratio mito volume / cytoplasm volume
            Rm_cell = Vmito # Volume ratio mito volume / cell volume
            Rc_cell = Vcyto # Volume ratio of cytosol / cell volume
        elif ExpType == 3: ## Diabetes
            W_c = 0.77 # cytoplasm water fraction - LeFurgey et al. 1991
            Vmito = 0.25*1.15 # (ml mito / ml cell)
            Vcyto = 0.678 # (ml cyto / ml cell)
            Rm_cyto = Vmito / Vcyto # Volume ratio mito volume / cytoplasm volume
            Rm_cell = Vmito # Volume ratio mito volume / cell volume
            Rc_cell = Vcyto # Volume ratio of cytosol / cell volume
        elif ExpType == 2: # Bose et al.'s experiments (JBC, 2003)
            W_c = 80 * 10 # Buffer water space
        else:
            W_c = 80
    elif tubule == "mTAL":
        if ExpType == 1:
            W_c = 0.77 # cytoplasm water fraction - LeFurgey et al. 1991
            Vmito = 0.33 # (ml mito / ml cell)
            Vcyto = 0.66 # (ml cyto / ml cell)
            Rm_cyto = Vmito / Vcyto # Volume ratio mito volume / cytoplasm volume
            Rm_cell = Vmito # Volume ratio mito volume / cell volume
            Rc_cell = Vcyto # Volume ratio of cytosol / cell volume
        elif ExpType == 2: # Bose et al.'s experiments (JBC, 2003)
            W_c = 80 * 10 # Buffer water space
        elif ExpType == 3: # Diabetes
            W_c = 0.77 # cytoplasm water fraction - LeFurgey et al. 1991
            Vmito = 0.33*1.15 # (ml mito / ml cell)
            Vcyto = 0.61 # (ml cyto / ml cell)
            Rm_cyto = Vmito / Vcyto # Volume ratio mito volume / cytoplasm volume
            Rm_cell = Vmito # Volume ratio mito volume / cell volume
            Rc_cell = Vcyto # Volume ratio of cytosol / cell volume
        else:
            W_c = 80

    ## Computing time derivatives of state variables

    f = np.zeros(pc.N_states)

    ## (i) Matrix species and dPsi
    f[pcIS.idPsi] = 1.48e5*( 4*Jsol1[pcFl.J_C1] + 2*Jsol1[pcFl.J_C3] +
        4*Jsol1[pcFl.J_C4] - pcPC.n_A*Jsol1[pcFl.J_F1] - Jsol1[pcFl.J_ANT]
        - Jsol1[pcFl.J_Hle] - Jsol1[pcFl.J_Kle] + Jsol1[pcFl.J_ASP_GLU])/\
        capacitanceWeight
    f[pcIS.iATP_x] = (Jsol1[pcFl.J_ndk] + Jsol1[pcFl.J_F1] - Jsol1[pcFl.J_ANT])/pcPC.W_x
    f[pcIS.iADP_x]  = (-Jsol1[pcFl.J_ndk] - Jsol1[pcFl.J_F1] + Jsol1[pcFl.J_ANT])/pcPC.W_x
    f[pcIS.iAMP_x]  = 0./pcPC.W_x
    f[pcIS.iGTP_x]  = (+Jsol1[pcFl.J_scoas] - Jsol1[pcFl.J_ndk])/pcPC.W_x
    f[pcIS.iGDP_x]  = (-Jsol1[pcFl.J_scoas] + Jsol1[pcFl.J_ndk])/pcPC.W_x
    f[pcIS.iPi_x]   = (-Jsol1[pcFl.J_scoas] - Jsol1[pcFl.J_F1] + Jsol1[pcFl.J_Pi1]
        - Jsol1[pcFl.J_MAL_PI])/pcPC.W_x
    f[pcIS.iNADH_x] = (+Jsol1[pcFl.J_pdh] + Jsol1[pcFl.J_isod] + Jsol1[pcFl.J_akgd] +
        Jsol1[pcFl.J_mdh] - Jsol1[pcFl.J_C1] + Jsol1[pcFl.J_gdh])/pcPC.W_x
    f[pcIS.iQH2_x]  = (+Jsol1[pcFl.J_sdh] + Jsol1[pcFl.J_C1] - Jsol1[pcFl.J_C3])/pcPC.W_x
    f[pcIS.iPYR_x]  = (-Jsol1[pcFl.J_pdh] + Jsol1[pcFl.J_PYR_H])/pcPC.W_x
    f[pcIS.iACCOA_x] = (-Jsol1[pcFl.J_cits] + Jsol1[pcFl.J_pdh])/pcPC.W_x
    f[pcIS.iCIT_x]  = (+Jsol1[pcFl.J_cits] - Jsol1[pcFl.J_acon] + Jsol1[pcFl.J_CIT_MAL])/ \
        pcPC.W_x
    f[pcIS.iICIT_x] = (+Jsol1[pcFl.J_acon] - Jsol1[pcFl.J_isod])/pcPC.W_x
    f[pcIS.iAKG_x]  = (+Jsol1[pcFl.J_isod] - Jsol1[pcFl.J_akgd] - Jsol1[pcFl.J_got] +
        Jsol1[pcFl.J_AKG_MAL] + Jsol1[pcFl.J_gdh])/pcPC.W_x
    f[pcIS.iSCOA_x] = (+Jsol1[pcFl.J_akgd] - Jsol1[pcFl.J_scoas])/pcPC.W_x
    f[pcIS.iCOASH_x]  = (-Jsol1[pcFl.J_pdh] - Jsol1[pcFl.J_akgd] + Jsol1[pcFl.J_scoas]
        + Jsol1[pcFl.J_cits])/pcPC.W_x
    f[pcIS.iSUC_x]  = (+Jsol1[pcFl.J_scoas] - Jsol1[pcFl.J_sdh] +
        Jsol1[pcFl.J_SUC_MAL])/pcPC.W_x
    f[pcIS.iFUM_x]  = (+Jsol1[pcFl.J_sdh] - Jsol1[pcFl.J_fum])/pcPC.W_x
    f[pcIS.iMAL_x]  = (+Jsol1[pcFl.J_fum] - Jsol1[pcFl.J_mdh] + Jsol1[pcFl.J_MAL_PI] -
        Jsol1[pcFl.J_AKG_MAL] - Jsol1[pcFl.J_CIT_MAL] - Jsol1[pcFl.J_SUC_MAL])/pcPC.W_x
    f[pcIS.iOAA_x]  = (-Jsol1[pcFl.J_cits] + Jsol1[pcFl.J_mdh] + Jsol1[pcFl.J_got])/pcPC.W_x

    f[pcIS.iGLU_x]  = (+Jsol1[pcFl.J_got] + Jsol1[pcFl.J_GLU_H]
        - Jsol1[pcFl.J_ASP_GLU] - Jsol1[pcFl.J_gdh])/pcPC.W_x
    f[pcIS.iASP_x]  = (-Jsol1[pcFl.J_got] + Jsol1[pcFl.J_ASP_GLU])/pcPC.W_x

    ## (ii) IM space species
    f[pcIS.iCred_i] = (+2*Jsol1[pcFl.J_C3] - 2*Jsol1[pcFl.J_C4])/pcPC.W_i
    f[pcIS.iATP_i] = (Jsol1[pcFl.J_ATP] + Jsol1[pcFl.J_ANT] + Jsol1[pcFl.J_AKi])/pcPC.W_i
    f[pcIS.iADP_i] = (Jsol1[pcFl.J_ADP] - Jsol1[pcFl.J_ANT] - 2*Jsol1[pcFl.J_AKi])/pcPC.W_i
    f[pcIS.iAMP_i] = (Jsol1[pcFl.J_AMP] + Jsol1[pcFl.J_AKi])/pcPC.W_i
    f[pcIS.iPi_i]  = (-Jsol1[pcFl.J_Pi1] + Jsol1[pcFl.J_Pi2] +
        Jsol1[pcFl.J_MAL_PI])/pcPC.W_i
    f[pcIS.iPYR_i] = (-Jsol1[pcFl.J_PYR_H] + Jsol1[pcFl.J_PYRt])/pcPC.W_i
    f[pcIS.iCIT_i] = (-Jsol1[pcFl.J_CIT_MAL] + Jsol1[pcFl.J_CITt])/pcPC.W_i
    f[pcIS.iICIT_i] = (Jsol1[pcFl.J_ICITt])/pcPC.W_i
    f[pcIS.iAKG_i] = (-Jsol1[pcFl.J_AKG_MAL] + Jsol1[pcFl.J_AKGt])/pcPC.W_i
    f[pcIS.iSUC_i] = (+ Jsol1[pcFl.J_SUCt] - Jsol1[pcFl.J_SUC_MAL])/pcPC.W_i
    f[pcIS.iFUM_i] = (Jsol1[pcFl.J_FUMt])/pcPC.W_i
    f[pcIS.iMAL_i] = (-Jsol1[pcFl.J_MAL_PI] + Jsol1[pcFl.J_MALt] +
        Jsol1[pcFl.J_AKG_MAL] + Jsol1[pcFl.J_CIT_MAL] +
        Jsol1[pcFl.J_SUC_MAL])/pcPC.W_i
    f[pcIS.iGLU_i] = (-Jsol1[pcFl.J_GLU_H] + Jsol1[pcFl.J_ASP_GLU] +
        Jsol1[pcFl.J_GLUt])/pcPC.W_i
    f[pcIS.iASP_i] = (-Jsol1[pcFl.J_ASP_GLU] + Jsol1[pcFl.J_ASPt])/pcPC.W_i

    ## (iii) Buffer species
    #  d[ADP]/dt, d[ATP]/dt and d[PI]/dt in cytoplasm are assumed to be
    #  consumed by Hexokinase reaction
    #  W_c   [=] buffer water/mito || cyto water/cyto
    #  J_ATP [=] M/s/l mito
    #  J_HK  [=] M/s/l buffer
    #  J_AtC, J_CKe [=] M/s/l cyto
    #  f [=] M/s/l cyto water
    if (ExpType == 2) or (ExpType == 5) or (ExpType == 6):
        f[pcIS.iPYR_c] = (-Jsol1[pcFl.J_PYRt])/W_c
        f[pcIS.iCIT_c] = (-Jsol1[pcFl.J_CITt])/W_c
        f[pcIS.iICIT_c] = (-Jsol1[pcFl.J_ICITt])/W_c
        f[pcIS.iAKG_c] = (-Jsol1[pcFl.J_AKGt])/W_c
        f[pcIS.iSUC_c] = (-Jsol1[pcFl.J_SUCt])/W_c
        f[pcIS.iFUM_c] = (-Jsol1[pcFl.J_FUMt])/W_c
        f[pcIS.iMAL_c] = (-Jsol1[pcFl.J_MALt])/W_c
        f[pcIS.iGLU_c] = (-Jsol1[pcFl.J_GLUt])/W_c
        f[pcIS.iASP_c] = (-Jsol1[pcFl.J_ASPt])/W_c
        f[pcIS.iPi_c]  = (-Jsol1[pcFl.J_Pi2])/W_c
        f[pcIS.iATP_c] = -Jsol1[pcFl.J_HK]/(W_c/(1+W_c))- \
            Jsol1[pcFl.J_ATP]/W_c # SET TO ZERO IN P/O MEASUREMENTS
        f[pcIS.iADP_c] = +Jsol1[pcFl.J_HK]/(W_c/(1+W_c))- \
            Jsol1[pcFl.J_ADP]/W_c # SET TO ZERO IN P/O MEASUREMENTS
        f[pcIS.iGLC_c] = -Jsol1[pcFl.J_HK]/(W_c/(1+W_c))
        f[pcIS.iG6P_c] = +Jsol1[pcFl.J_HK]/(W_c/(1+W_c))
    else:
        f[pcIS.iPYR_c] = (-Rm_cyto*Jsol1[pcFl.J_PYRt])/W_c
        f[pcIS.iCIT_c] = (-Rm_cyto*Jsol1[pcFl.J_CITt])/W_c
        f[pcIS.iICIT_c] = (-Rm_cyto*Jsol1[pcFl.J_ICITt])/W_c
        f[pcIS.iAKG_c] = (-Rm_cyto*Jsol1[pcFl.J_AKGt])/W_c
        f[pcIS.iSUC_c] = (-Rm_cyto*Jsol1[pcFl.J_SUCt])/W_c
        f[pcIS.iFUM_c] = (-Rm_cyto*Jsol1[pcFl.J_FUMt])/W_c
        f[pcIS.iMAL_c] = (-Rm_cyto*Jsol1[pcFl.J_MALt])/W_c
        f[pcIS.iGLU_c] = (-Rm_cyto*Jsol1[pcFl.J_GLUt])/W_c
        f[pcIS.iASP_c] = (-Rm_cyto*Jsol1[pcFl.J_ASPt])/W_c
        f[pcIS.iPi_c]  = (-Rm_cyto*Jsol1[pcFl.J_Pi2]+J_AtC)/W_c
        f[pcIS.iATP_c] = (-Rm_cyto*Jsol1[pcFl.J_ATP]-J_AtC+
            Jsol1[pcFl.J_CKe]+Jsol1[pcFl.J_AKe])/W_c
        f[pcIS.iADP_c] = (-Rm_cyto*Jsol1[pcFl.J_ADP]+J_AtC-
            Jsol1[pcFl.J_CKe]-2*Jsol1[pcFl.J_AKe])/W_c
        f[pcIS.iPCr_c] = (-Jsol1[pcFl.J_CKe])/W_c
        f[pcIS.iAMP_c] = (-Rm_cyto*Jsol1[pcFl.J_AMP]+Jsol1[pcFl.J_AKe])/W_c

    if (StateType == 0):
        f[pcIS.iPYR_c] = 0.
        f[pcIS.iMAL_c] = 0.
        f[pcIS.iCIT_c] = 0.
        f[pcIS.iICIT_c] = 0.
        f[pcIS.iAKG_c] = 0.
        f[pcIS.iSUC_c] = 0.
        f[pcIS.iFUM_c] = 0.

    if (ExpType == 1) or (ExpType == 3) or (ExpType == 4):
        f[pcIS.iPYR_c] = 0.
    elif (ExpType == 5): # for uncoupling
        f[pcIS.iPYR_c] = 0.
        f[pcIS.iMAL_c] = 0.
        f[pcIS.iCIT_c] = 0.
        f[pcIS.iICIT_c] = 0.
        f[pcIS.iAKG_c] = 0.
        f[pcIS.iSUC_c] = 0.
        f[pcIS.iFUM_c] = 0.

    # %% Computing dH/dt, dMg/dt, and dK/dt
    # % Calculate dH/dt, dMg/dt, dK/dt according to Dan's strategy (available
    # % in the chapter "Biochemical Reaction Networks" in Beard and Qian's book)
    #
    # % assume [H+], [Mg2+], [K+] constant in IM and cytoplasm/buffer space
    # % all the below calculation if for dH_x/dt, dMg_x/dt, dK_x/dt
    #
    # % take concentration of reactants located in the mito. matrix
    # % from state variables

    xx = np.zeros(pc.N_reactant)
    xx[pcIR.iH] = x[pcIS.iH_x]
    xx[pcIR.iATP] = x[pcIS.iATP_x]
    xx[pcIR.iADP] = x[pcIS.iADP_x]
    xx[pcIR.iAMP] = x[pcIS.iAMP_x]
    xx[pcIR.iGTP] = x[pcIS.iGTP_x]
    xx[pcIR.iGDP] = x[pcIS.iGDP_x]
    xx[pcIR.iPi] = x[pcIS.iPi_x]
    xx[pcIR.iNADH] = x[pcIS.iNADH_x]
    xx[pcIR.iNAD] = pcPC.NADtot - x[pcIS.iNADH_x]
    xx[pcIR.iQH2] = x[pcIS.iQH2_x]
    xx[pcIR.iCOQ] = pcPC.Qtot - x[pcIS.iQH2_x]
    xx[pcIR.iPYR] = x[pcIS.iPYR_x]
    xx[pcIR.iOAA] = x[pcIS.iOAA_x]
    xx[pcIR.iACCOA] = x[pcIS.iACCOA_x]
    xx[pcIR.iCIT] = x[pcIS.iCIT_x]
    xx[pcIR.iICIT] = x[pcIS.iICIT_x]
    xx[pcIR.iAKG] = x[pcIS.iAKG_x]
    xx[pcIR.iSCOA] = x[pcIS.iSCOA_x]
    xx[pcIR.iCOASH] = x[pcIS.iCOASH_x]
    xx[pcIR.iSUC] = x[pcIS.iSUC_x]
    xx[pcIR.iFUM] = x[pcIS.iFUM_x]
    xx[pcIR.iMAL] = x[pcIS.iMAL_x]
    xx[pcIR.iGLU] = x[pcIS.iGLU_x]
    xx[pcIR.iASP] = x[pcIS.iASP_x]
    xx[pcIR.iK] = x[pcIS.iK_x]
    xx[pcIR.iMg] = x[pcIS.iMg_x]
    xx[pcIR.iO2] = x[pcIS.iO2_x]
    xx[pcIR.iFADH2] = pcPC.FADtot/2
    xx[pcIR.iFAD] = pcPC.FADtot - xx[pcIR.iFADH2]
    xx[pcIR.iCO2tot] = x[pcIS.iCO2tot_x]

    dxxdt = np.zeros(pc.N_reactant)
    # dxxdt[pcIR.iH] = f[pcIS.iH_x]
    dxxdt[pcIR.iATP] = f[pcIS.iATP_x]
    dxxdt[pcIR.iADP] = f[pcIS.iADP_x]
    dxxdt[pcIR.iAMP] = f[pcIS.iAMP_x]
    dxxdt[pcIR.iGTP] = f[pcIS.iGTP_x]
    dxxdt[pcIR.iGDP] = f[pcIS.iGDP_x]
    dxxdt[pcIR.iPi] = f[pcIS.iPi_x]
    dxxdt[pcIR.iNADH] = f[pcIS.iNADH_x]
    dxxdt[pcIR.iNAD] = -f[pcIS.iNADH_x]
    dxxdt[pcIR.iQH2] = f[pcIS.iQH2_x]
    dxxdt[pcIR.iCOQ] = - f[pcIS.iQH2_x]
    dxxdt[pcIR.iPYR] = f[pcIS.iPYR_x]
    dxxdt[pcIR.iOAA] = f[pcIS.iOAA_x]
    dxxdt[pcIR.iACCOA] = f[pcIS.iACCOA_x]
    dxxdt[pcIR.iCIT] = f[pcIS.iCIT_x]
    dxxdt[pcIR.iICIT] = f[pcIS.iICIT_x]
    dxxdt[pcIR.iAKG] = f[pcIS.iAKG_x]
    dxxdt[pcIR.iSCOA] = f[pcIS.iSCOA_x]
    dxxdt[pcIR.iCOASH] = f[pcIS.iCOASH_x]
    dxxdt[pcIR.iSUC] = f[pcIS.iSUC_x]
    dxxdt[pcIR.iFUM] = f[pcIS.iFUM_x]
    dxxdt[pcIR.iMAL] = f[pcIS.iMAL_x]
    dxxdt[pcIR.iGLU] = f[pcIS.iGLU_x]
    dxxdt[pcIR.iASP] = f[pcIS.iASP_x]
    # dxxdt[pcIR.iK]   = f[pcIS.iK_x]
    # dxxdt[pcIR.iMg]  = f[pcIS.iMg_x]
    dxxdt[pcIR.iCred] = f[pcIS.iCred_i]
    dxxdt[pcIR.iCox] = -f[pcIS.iCred_i]
    dxxdt[pcIR.iO2]  = f[pcIS.iO2_x]
    # dxxdt[pcIR.iH2O] = 0
    # dxxdt[pcIR.iFADH2] = 0
    dxxdt[pcIR.iFAD] = -dxxdt[pcIR.iFADH2]
    dxxdt[pcIR.iCO2tot] = f[pcIS.iCO2tot_x]

    # Necessary concentrations
    H_x = x[pcIS.iH_x]
    K_x = x[pcIS.iK_x]
    Mg_x = x[pcIS.iMg_x]

    ## K values for reference species
    ## Repeated code from fluxes.py
    Kh = [float('inf')] * pc.N_reactant
    Km = [float('inf')] * pc.N_reactant
    Kk = [float('inf')] * pc.N_reactant
    Kh[pcIR.iATP] = 10 ** (-6.59)
    Km[pcIR.iATP] = 10 ** (-3.82)
    Kk[pcIR.iATP] = 10 ** (-1.013)
    Kh[pcIR.iADP] = 10 ** (-6.42)
    Km[pcIR.iADP] = 10 ** (-2.79)
    Kk[pcIR.iADP] = 10 ** (-0.882)
    Kh[pcIR.iAMP] = 10 ** (-6.22)
    Km[pcIR.iAMP] = 10 ** (-1.86)
    Kk[pcIR.iAMP] = 10 ** (-0.6215)
    Kh[pcIR.iGTP] = Kh[pcIR.iATP]
    Km[pcIR.iGTP] = Km[pcIR.iATP]
    Kk[pcIR.iGTP] = Kk[pcIR.iATP]
    Kh[pcIR.iGDP] = Kh[pcIR.iADP]
    Km[pcIR.iGDP] = Km[pcIR.iADP]
    Kk[pcIR.iGDP] = Kk[pcIR.iADP]
    Kh[pcIR.iPi] = 10 ** (-6.71)
    Km[pcIR.iPi] = 10 ** (-1.69)
    Kk[pcIR.iPi] = 10 ** (+0.0074)
    Kh[pcIR.iCOASH] = 10 ** (-8.13)
    Km[pcIR.iOAA] = 10 ** (-0.8629)
    Kh[pcIR.iCIT] = 10 ** (-5.63)
    Km[pcIR.iCIT] = 10 ** (-3.37)
    Kk[pcIR.iCIT] = 10 ** (-0.339)
    Kh[pcIR.iICIT] = 10 ** (-5.64)
    Km[pcIR.iICIT] = 10 ** (-2.46)
    Kh[pcIR.iSCOA] = 10 ** (-3.96)
    Kh[pcIR.iSUC] = 10 ** (-5.13)
    Km[pcIR.iSUC] = 10 ** (-1.17)
    Kk[pcIR.iSUC] = 10 ** (-0.3525)
    Kh[pcIR.iFUM] = 10 ** (-4.10)
    Kh[pcIR.iMAL] = 10 ** (-4.75)
    Km[pcIR.iMAL] = 10 ** (-1.55)
    Kk[pcIR.iMAL] = 10 ** (+0.107)
    Kh[pcIR.iCO2tot] = 10 ** (-9.82)
    Km[pcIR.iPYR] = 10 ** (-1.02)
    Kh[pcIR.iG6P] = 10 ** (-5.91)
    # Kh(iGLU) = 10**(-4.25) # from Nelson & Cox, "Lehninger's Princinples of Biochemistry"
    # , p78
    Kh[pcIR.iGLU] = 10 ** (-4.06)  # 37 C, I = 0.15
    Km[pcIR.iGLU] = 10 ** (-1.82)
    # Kh(iASP) = 10 ^ (-3.65) # from Nelson & Cox, "Lehninger's Princinples of Biochemistry"
    # , p78
    Kh[pcIR.iASP] = 10 ** (-3.65)  # 37 C, I = 0.15
    Km[pcIR.iASP] = 10 ** (-2.32)

    ## # compute binding polynomials for reactants, also repeated code from fluxes.py
    P_x = [1.]*pc.N_reactant
    P_x[pcIR.iATP] = 1 + H_x/Kh[pcIR.iATP] + Mg_x/Km[pcIR.iATP] + K_x/Kk[pcIR.iATP]
    P_x[pcIR.iADP] = 1 + H_x/Kh[pcIR.iADP] + Mg_x/Km[pcIR.iADP] + K_x/Kk[pcIR.iADP]
    P_x[pcIR.iAMP] = 1 + H_x/Kh[pcIR.iAMP] + Mg_x/Km[pcIR.iAMP] + K_x/Kk[pcIR.iAMP]
    P_x[pcIR.iGTP] = P_x[pcIR.iATP]
    P_x[pcIR.iGDP] = P_x[pcIR.iADP]
    P_x[pcIR.iPi] = 1 + H_x/Kh[pcIR.iPi] + Mg_x/Km[pcIR.iPi] + K_x/Kk[pcIR.iPi]
    P_x[pcIR.iCOASH] = 1 + H_x/Kh[pcIR.iCOASH]
    P_x[pcIR.iOAA] = 1 + Mg_x/Km[pcIR.iOAA]
    P_x[pcIR.iCIT] = 1 + H_x/Kh[pcIR.iCIT] + Mg_x/Km[pcIR.iCIT] + K_x/Kk[pcIR.iCIT]
    P_x[pcIR.iICIT] = 1 + H_x/Kh[pcIR.iICIT] + Mg_x/Km[pcIR.iICIT]
    P_x[pcIR.iSCOA] = 1 + H_x/Kh[pcIR.iSCOA]
    P_x[pcIR.iSUC] = 1 + H_x/Kh[pcIR.iSUC] + Mg_x/Km[pcIR.iSUC] + K_x/Kk[pcIR.iSUC]
    P_x[pcIR.iFUM] = 1 + H_x/Kh[pcIR.iFUM]
    P_x[pcIR.iMAL] = 1 + H_x/Kh[pcIR.iMAL] + Mg_x/Km[pcIR.iMAL] + K_x/Kk[pcIR.iMAL]
    P_x[pcIR.iCO2tot] = 1 + H_x/Kh[pcIR.iCO2tot]
    P_x[pcIR.iPYR] = 1 + Mg_x/Km[pcIR.iPYR]
    P_x[pcIR.iG6P] = 1 + H_x/Kh[pcIR.iG6P]
    P_x[pcIR.iGLU] = 1 + H_x/Kh[pcIR.iGLU] + Mg_x/Km[pcIR.iGLU]
    P_x[pcIR.iASP] = 1 + H_x/Kh[pcIR.iASP] + Mg_x/Km[pcIR.iASP]

    P_x = np.array(P_x)
    H_x = np.array(H_x)
    Mg_x = np.array(Mg_x)
    K_x = np.array(K_x)
    Kh = np.array(Kh)
    Km = np.array(Km)
    Kk = np.array(Kk)

    # Partial Derivatives
    pHBpM = - np.sum( (H_x*xx/Kh)/(Km*np.power(P_x,2)) )
    pHBpK = - np.sum( (H_x*xx/Kh)/(Kk*np.power(P_x,2)) )
    pHBpH = np.sum( (1+Mg_x/Km+K_x/Kk)*xx/(Kh*np.power(P_x,2)) )
    pMBpH = - np.sum( (Mg_x*xx/Km)/(Kh*np.power(P_x,2)) )
    pMBpK = - np.sum( (Mg_x*xx/Km)/(Kk*np.power(P_x,2)) )
    pMBpM = np.sum( (1+H_x/Kh+K_x/Kk)*xx/(Km*np.power(P_x,2)) )
    pKBpH = - np.sum( (K_x*xx/Kk)/(Kh*np.power(P_x,2)) )
    pKBpM = - np.sum( (K_x*xx/Kk)/(Km*np.power(P_x,2)) )
    pKBpK = np.sum( (1+H_x/Kh+Mg_x/Km)*xx/(Kk*np.power(P_x,2)) )

    # PHI's:
    Phi_H = - np.sum( H_x*dxxdt/(Kh*P_x) ) + 1*(-1*Jsol1[pcFl.J_pdh] +
        2*Jsol1[pcFl.J_cits] + (-1)*Jsol1[pcFl.J_akgd] + Jsol1[pcFl.J_scoas]
        + Jsol1[pcFl.J_mdh]+ 1*(Jsol1[pcFl.J_PYR_H] + Jsol1[pcFl.J_GLU_H] +
        Jsol1[pcFl.J_CIT_MAL] - Jsol1[pcFl.J_ASP_GLU])-(4+1)*Jsol1[pcFl.J_C1] -
        (4-2)*Jsol1[pcFl.J_C3] - (2+2)*Jsol1[pcFl.J_C4] + (pcPC.n_A-1)*Jsol1[pcFl.J_F1]
        + 2*Jsol1[pcFl.J_Pi1] + 1*Jsol1[pcFl.J_Hle] - 1*Jsol1[pcFl.J_KH]
        + Jsol1[pcFl.J_gdh])/pcPC.W_x
    Phi_M = -np.sum( Mg_x*dxxdt/(Km*P_x) )
    Phi_K = - np.sum( K_x*dxxdt/(Kk*P_x) ) + 1*Jsol1[pcFl.J_Kle]/pcPC.W_x \
        + 1*Jsol1[pcFl.J_KH]/pcPC.W_x

    # alpha's:
    aH = 1 + pHBpH
    aM = 1 + pMBpM
    aK = 1 + pKBpK

    # add additional buffer for [H+]
    BX = 0.02 # M
    K_BX = 1e-7 # M
    aH = 1 + pHBpH + BX/K_BX/(1+H_x/K_BX)**2

    #  Denominator:
    D = aH*pKBpM*pMBpK + aK*pHBpM*pMBpH + aM*pHBpK*pKBpH - \
        aM*aK*aH - pHBpK*pKBpM*pMBpH - pHBpM*pMBpK*pKBpH

    # Derivatives for H, Mg, K:
    dH_xdt = ((pKBpM * pMBpK - aM * aK) * Phi_H +
        (aK * pHBpM - pHBpK * pKBpM) * Phi_M +
        (aM * pHBpK - pHBpM * pMBpK) * Phi_K) / D

    dMg_xdt = ((aK * pMBpH - pKBpH * pMBpK) * Phi_H +
        (pKBpH * pHBpK - aH * aK) * Phi_M +
        (aH * pMBpK - pHBpK * pMBpH) * Phi_K) / D

    dK_xdt = ((aM * pKBpH - pKBpM * pMBpH) * Phi_H +
        (aH * pKBpM - pKBpH * pHBpM) * Phi_M +
        (pMBpH * pHBpM - aH * aM) * Phi_K) / D

    f[pcIS.iH_x] = dH_xdt
    f[pcIS.iMg_x] = dMg_xdt
    f[pcIS.iK_x] = dK_xdt

    return f

def conservationEqsmTAL(x, J_AtC, ExpType = 1, StateType = 1,
                    w = [1., 1., 1., 1.], timeStart = 1,
                    capacitanceWeight = 1.,
                    potassiumW = 1.,
                    glyc = 0.0, DCA = 1., cana = 1.,
                    params = paramsmTAL):
    return conservationEqs1(x = x, J_AtC = J_AtC, ExpType = ExpType,
                            StateType = StateType,
                            w = w, timeStart = timeStart,
                            capacitanceWeight = capacitanceWeight,
                            potassiumW = potassiumW, tubule = "mTAL",
                            glyc = glyc, DCA = DCA, cana = cana,
                            params = paramsmTAL)