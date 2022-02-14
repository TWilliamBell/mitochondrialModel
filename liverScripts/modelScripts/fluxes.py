import numpy as np

from pc import * ## Contains the 'namespaces' pcITA, pcIR, pcIS, pcPC, pcITCA, and
## the list ics.

def fluxes(x, param, ExpType, w = [1., 1., 1., 1.], potassiumW = 1., DCA = 1.):
    """This function calculates all the Mitochondrial fluxes for the model.
    It is based on similar code copyrighted by Fan Wu and Daniel. A Beard,
    and subsequent modifications of that code by Edwards et al. 2019.
    The ouput is all the fluxes of the system.  The inputs are x - the state
    of the cell, param - the parameter values used, and ExpType - specifies whether
    the process is being done in a bath or in vivo (1 in vivo, 2 in buffer,
    5 uncoupling)."""

    # For anyone comparing with the original code, the indices of the relevant
    # concentrations have been moved to the pc namespace.

    # Vmax values for TCA fluxes from param
    Vmax = [0.]*11
    Vmax[0] = param[0]
    Vmax[1] = param[1] # 0.199; CS activity - Benard et al., 2016
    Vmax[2] = param[2] # Above comment from the original code.
    Vmax[3] = param[3]
    Vmax[4] = param[4]
    Vmax[5] = param[5]
    Vmax[6] = param[6]
    Vmax[7] = param[7]
    Vmax[8] = param[8]
    Vmax[9] = param[9]
    Vmax[10] = param[10]

    # Activities of transporters, parameters 15, 18, 19, 21-23, 25-29 not used
    x_PYR_H = param[11]
    x_GLU_H = param[12]
    x_CIT_MAL = param[13]
    x_AKG_MAL = param[14]
    x_MAL_PI  = param[16]
    x_ASP_GLU = param[17]
    x_SUC_MAL = param[20]

    # Hexokinase activity is not included
    x_HK = 0

    # Parameters for oxidative phosphorylation
    x_C1 = param[30] * 1.2
    x_C3 = param[31] * 0.525
    x_C4 = param[32] * 0.40
    x_F1 = param[33]
    x_ANT = param[34] # Best - fit 0.85 before
    x_Pi1 = param[35]
    k_PiH = param[36]
    x_KH = param[37]
    x_Hle = param[38] # *0.61; % Best - fit 0.61 prior to leak proton change
    k_Pi1 = param[39]
    k_Pi2 = param[40]

    P_Hle = 345.00 * param[38] / 302.0928 * 1 # 345 used to be 750 in absence
    # of Potassium leak
    P_Kle = 345.00 * 1.e-006 * 1 * potassiumW
    x_Kle = x_Hle * 0 # 0.5e-006; #2.0e-006;

    # Medulla vs. cortex simulations
    # x_C1 = param[30] * 0.75 * 2.0
    # x_C3 = param[31] * 0.50 * 2.0 # 0.516 From Benard et al., 2006
    # x_C4 = param[32] * 0.25 * 2.0 # 0.264 From Benard et al., 2006
    # x_F1 = param[33] * 2.0
    # x_Pi1 = param[35] * 2.0

    # Permeability coefficients for TCA intermediates
    p_TI = 85 # assumed to be equal to x_A for nucleotides(micron sec ^ -1)
    x_PYRt = p_TI
    x_GLUt = p_TI
    x_ASPt = p_TI
    x_CITt = p_TI
    x_ICITt = p_TI
    x_AKGt = p_TI
    x_FUMt = 0 * p_TI
    x_SUCt = p_TI
    x_MALt = p_TI

    ## Pooled concentrations and OXPHOS-related quantities moved to pcPC namespace
    ## in pc.py because they're used in more than one script.

    # Outer membrane transport parameters
    x_A = 85 # micron sec ^ {-1}
    x_Pi2 = 327 # micron sec ^ {-1}
    gamma = 5.99 # mito membrane area per cell volume micron ^ {-1}

    # Load values of state variables

    MinCon = 1e-32

    # Matrix species and dPsi

    dPsi = x[pcIS.idPsi]
    ATP_x = x[pcIS.iATP_x]
    ADP_x = x[pcIS.iADP_x]
    AMP_x = x[pcIS.iAMP_x]
    GTP_x = x[pcIS.iGTP_x]
    GDP_x = x[pcIS.iGDP_x]
    Pi_x = x[pcIS.iPi_x]
    NADH_x = x[pcIS.iNADH_x]
    QH2_x = x[pcIS.iQH2_x]
    PYR_x = x[pcIS.iPYR_x]
    OAA_x = x[pcIS.iOAA_x]
    ACCOA_x = x[pcIS.iACCOA_x]
    CIT_x = x[pcIS.iCIT_x]
    ICIT_x = x[pcIS.iICIT_x]
    AKG_x = x[pcIS.iAKG_x]
    SCOA_x = x[pcIS.iSCOA_x]
    COASH_x = x[pcIS.iCOASH_x]
    SUC_x = x[pcIS.iSUC_x]
    FUM_x = x[pcIS.iFUM_x]
    MAL_x = x[pcIS.iMAL_x]
    GLU_x = x[pcIS.iGLU_x]
    ASP_x = x[pcIS.iASP_x]
    H_x = x[pcIS.iH_x]
    K_x = x[pcIS.iK_x]
    Mg_x = x[pcIS.iMg_x]
    O2 = x[pcIS.iO2_x]
    CO2tot = x[pcIS.iCO2tot_x]

    #(ii) IM space species
    Cred_i = max(0, x[pcIS.iCred_i])
    ATP_i = x[pcIS.iATP_i]
    ADP_i = x[pcIS.iADP_i]
    AMP_i = x[pcIS.iAMP_i]
    Pi_i = x[pcIS.iPi_i]
    PYR_i = x[pcIS.iPYR_i]
    CIT_i = x[pcIS.iCIT_i]
    AKG_i = x[pcIS.iAKG_i]
    SUC_i = x[pcIS.iSUC_i]
    MAL_i = x[pcIS.iMAL_i]
    GLU_i = x[pcIS.iGLU_i]
    ASP_i = x[pcIS.iASP_i]
    H_i = x[pcIS.iH_i]
    Mg_i = x[pcIS.iMg_i]
    K_i = x[pcIS.iK_i]
    FUM_i = x[pcIS.iFUM_i]
    ICIT_i = x[pcIS.iICIT_i]

    # (iii) Cytoplasmic species
    ATP_c = x[pcIS.iATP_c]
    ADP_c = x[pcIS.iADP_c]
    Pi_c = x[pcIS.iPi_c]
    PYR_c = x[pcIS.iPYR_c]
    CIT_c = x[pcIS.iCIT_c]
    AKG_c = x[pcIS.iAKG_c]
    SUC_c = x[pcIS.iSUC_c]
    MAL_c = x[pcIS.iMAL_c]
    GLU_c = x[pcIS.iGLU_c]
    ASP_c = x[pcIS.iASP_c]
    H_c = x[pcIS.iH_c]
    Mg_c = x[pcIS.iMg_c]
    K_c = x[pcIS.iK_c]
    FUM_c = x[pcIS.iFUM_c]
    ICIT_c = x[pcIS.iICIT_c]
    GLC_c = x[pcIS.iGLC_c]
    G6P_c = x[pcIS.iG6P_c]
    AMP_c = x[pcIS.iAMP_c]

    # (iv) Other concentrations computed from the state variables:
    NAD_x  = pcPC.NADtot - NADH_x
    COQ_x  = pcPC.Qtot - QH2_x
    Cox_i  = pcPC.Ctot - Cred_i

    # (v) set the H +, Mg2 +, and K + are permeable for outer mito membrane
    H_i = H_c
    Mg_i = Mg_c
    K_i = K_c

    # Loading thermodynamic data (deltG, pK, etc.)
    # T = 298.15K (25 C)    I = 0.17 M
    # Standard Gibbs free energy of formation of reference species, (kJ/mol)
    # without temperature correction on dGf

    dGf1 = [0.]*pc.N_reactant
    dGf1[pcIR.iH2O] = -235.74 # H2O
    dGf1[pcIR.iO2] = 16.40 # O2(aq)
    dGf1[pcIR.iNADH] = 39.31 # NADH
    dGf1[pcIR.iNAD] = 18.10 # NAD +
    dGf1[pcIR.iQH2] = -23.30 # QH2
    dGf1[pcIR.iCOQ] = 65.17 # Q
    dGf1[pcIR.iATP] = -2771.00 #  ATP4 -
    dGf1[pcIR.iADP] = -1903.96 # ADP3 -
    dGf1[pcIR.iAMP] = -1034.66 # AMP2 -
    dGf1[pcIR.iGTP] = dGf1[pcIR.iATP]
    dGf1[pcIR.iGDP] = dGf1[pcIR.iADP]
    dGf1[pcIR.iCred] = -27.41 # CytoC(red) 2 +
    dGf1[pcIR.iCox] = -6.52 # CytoC(ox) 3 +
    dGf1[pcIR.iPi] = -1098.27 # HPO42 -
    dGf1[pcIR.iPCr] = 0 # PCr2 -
    dGf1[pcIR.iCr] = -252.68 # HCr
    dGf1[pcIR.iFADH2] = -67.60 # FADH2 - enz
    dGf1[pcIR.iFAD] = 19.55 # FAD - enz
    dGf1[pcIR.iCOASH] = -0.72 # CoAS -
    dGf1[pcIR.iACCOA] = -178.19 # AcCoA
    dGf1[pcIR.iOAA] = -794.74 # OAA2 -
    dGf1[pcIR.iCIT] = -1165.59 # CIT3 -
    dGf1[pcIR.iICIT] = -1158.94 # ICIT3 -
    dGf1[pcIR.iAKG] = -793.41 # AKG2 -
    dGf1[pcIR.iSCOA] = -507.55 # SCoA -
    dGf1[pcIR.iSUC] = -690.44 # SUC2 -
    dGf1[pcIR.iFUM] = -603.32 # FUM2 -
    dGf1[pcIR.iMAL] = -842.66 # MAL2 -
    dGf1[pcIR.iASP] = -692.26 # ASP -
    dGf1[pcIR.iGLU] = -692.40 # GLU - (L - glutamate)
    dGf1[pcIR.iCO2tot] = -530.71 # CO2tot
    dGf1[pcIR.iPYR] = -470.82 # PYR2 -
    dGf1[pcIR.iGLC] = -907.21 # Glucose
    dGf1[pcIR.iG6P] = -1758.87 # Glucose - 6 - phosphate

    # K values for reference species
    Kh = [float('inf')]*pc.N_reactant
    Km = [float('inf')]*pc.N_reactant
    Kk = [float('inf')]*pc.N_reactant
    Kh[pcIR.iATP] = 10**(-6.59)
    Km[pcIR.iATP] = 10**(-3.82)
    Kk[pcIR.iATP] = 10**(-1.013)
    Kh[pcIR.iADP] = 10**(-6.42)
    Km[pcIR.iADP] = 10**(-2.79)
    Kk[pcIR.iADP] = 10**(-0.882)
    Kh[pcIR.iAMP] = 10**(-6.22)
    Km[pcIR.iAMP] = 10**(-1.86)
    Kk[pcIR.iAMP] = 10**(-0.6215)
    Kh[pcIR.iGTP] = Kh[pcIR.iATP]
    Km[pcIR.iGTP] = Km[pcIR.iATP]
    Kk[pcIR.iGTP] = Kk[pcIR.iATP]
    Kh[pcIR.iGDP] = Kh[pcIR.iADP]
    Km[pcIR.iGDP] = Km[pcIR.iADP]
    Kk[pcIR.iGDP] = Kk[pcIR.iADP]
    Kh[pcIR.iPi] = 10**(-6.71)
    Km[pcIR.iPi] = 10**(-1.69)
    Kk[pcIR.iPi] = 10**(+0.0074)
    Kh[pcIR.iCOASH] = 10**(-8.13)
    Km[pcIR.iOAA] = 10**(-0.8629)
    Kh[pcIR.iCIT] = 10**(-5.63)
    Km[pcIR.iCIT] = 10**(-3.37)
    Kk[pcIR.iCIT] = 10**(-0.339)
    Kh[pcIR.iICIT] = 10**(-5.64)
    Km[pcIR.iICIT] = 10**(-2.46)
    Kh[pcIR.iSCOA] = 10**(-3.96)
    Kh[pcIR.iSUC] = 10**(-5.13)
    Km[pcIR.iSUC] = 10**(-1.17)
    Kk[pcIR.iSUC] = 10**(-0.3525)
    Kh[pcIR.iFUM] = 10**(-4.10)
    Kh[pcIR.iMAL] = 10**(-4.75)
    Km[pcIR.iMAL] = 10**(-1.55)
    Kk[pcIR.iMAL] = 10**(+0.107)
    Kh[pcIR.iCO2tot] = 10**(-9.82)
    Km[pcIR.iPYR] = 10**(-1.02)
    Kh[pcIR.iG6P] = 10**(-5.91)
    # Kh(iGLU) = 10**(-4.25) # from Nelson & Cox, "Lehninger's Princinples of Biochemistry"
    # , p78
    Kh[pcIR.iGLU] = 10**(-4.06) # 37 C, I = 0.15
    Km[pcIR.iGLU] = 10**(-1.82)
    # Kh(iASP) = 10 ^ (-3.65) # from Nelson & Cox, "Lehninger's Princinples of Biochemistry"
    # , p78
    Kh[pcIR.iASP] = 10**(-3.65) # 37 C, I = 0.15
    Km[pcIR.iASP] = 10**(-2.32)

    # compute binding polynomials for reactants
    P_x = [1.]*pc.N_reactant
    P_c = [1.]*pc.N_reactant
    P_i = [1.]*pc.N_reactant
    P_x[pcIR.iATP] = 1 + H_x/Kh[pcIR.iATP] + Mg_x/Km[pcIR.iATP] + K_x/Kk[pcIR.iATP]
    P_c[pcIR.iATP] = 1 + H_c/Kh[pcIR.iATP] + Mg_c/Km[pcIR.iATP] + K_c/Kk[pcIR.iATP]
    P_i[pcIR.iATP] = 1 + H_i/Kh[pcIR.iATP] + Mg_i/Km[pcIR.iATP] + K_i/Kk[pcIR.iATP]
    P_x[pcIR.iADP] = 1 + H_x/Kh[pcIR.iADP] + Mg_x/Km[pcIR.iADP] + K_x/Kk[pcIR.iADP]
    P_c[pcIR.iADP] = 1 + H_c/Kh[pcIR.iADP] + Mg_c/Km[pcIR.iADP] + K_c/Kk[pcIR.iADP]
    P_i[pcIR.iADP] = 1 + H_i/Kh[pcIR.iADP] + Mg_i/Km[pcIR.iADP] + K_i/Kk[pcIR.iADP]
    P_x[pcIR.iAMP] = 1 + H_x/Kh[pcIR.iAMP] + Mg_x/Km[pcIR.iAMP] + K_x/Kk[pcIR.iAMP]
    P_c[pcIR.iAMP] = 1 + H_c/Kh[pcIR.iAMP] + Mg_c/Km[pcIR.iAMP] + K_c/Kk[pcIR.iAMP]
    P_i[pcIR.iAMP] = 1 + H_i/Kh[pcIR.iAMP] + Mg_i/Km[pcIR.iAMP] + K_i/Kk[pcIR.iAMP]
    P_x[pcIR.iGTP] = P_x[pcIR.iATP]
    P_c[pcIR.iGTP] = P_c[pcIR.iATP]
    P_i[pcIR.iGTP] = P_i[pcIR.iATP]
    P_x[pcIR.iGDP] = P_x[pcIR.iADP]
    P_c[pcIR.iGDP] = P_c[pcIR.iADP]
    P_i[pcIR.iGDP] = P_i[pcIR.iADP]
    P_x[pcIR.iPi] = 1 + H_x/Kh[pcIR.iPi] + Mg_x/Km[pcIR.iPi] + K_x/Kk[pcIR.iPi]
    P_c[pcIR.iPi] = 1 + H_c/Kh[pcIR.iPi] + Mg_c/Km[pcIR.iPi] + K_c/Kk[pcIR.iPi]
    P_i[pcIR.iPi] = 1 + H_i/Kh[pcIR.iPi] + Mg_i/Km[pcIR.iPi] + K_i/Kk[pcIR.iPi]
    P_x[pcIR.iCOASH] = 1 + H_x/Kh[pcIR.iCOASH]
    P_i[pcIR.iCOASH] = 1 + H_i/Kh[pcIR.iCOASH]
    P_c[pcIR.iCOASH] = 1 + H_c/Kh[pcIR.iCOASH]
    P_x[pcIR.iOAA] = 1 + Mg_x/Km[pcIR.iOAA]
    P_i[pcIR.iOAA] = 1 + Mg_i/Km[pcIR.iOAA]
    P_c[pcIR.iOAA] = 1 + Mg_c/Km[pcIR.iOAA]
    P_x[pcIR.iCIT] = 1 + H_x/Kh[pcIR.iCIT] + Mg_x/Km[pcIR.iCIT] + K_x/Kk[pcIR.iCIT]
    P_i[pcIR.iCIT] = 1 + H_i/Kh[pcIR.iCIT] + Mg_i/Km[pcIR.iCIT] + K_i/Kk[pcIR.iCIT]
    P_c[pcIR.iCIT] = 1 + H_c/Kh[pcIR.iCIT] + Mg_c/Km[pcIR.iCIT] + K_c/Kk[pcIR.iCIT]
    P_x[pcIR.iICIT] = 1 + H_x/Kh[pcIR.iICIT] + Mg_x/Km[pcIR.iICIT]
    P_i[pcIR.iICIT] = 1 + H_i/Kh[pcIR.iICIT] + Mg_i/Km[pcIR.iICIT]
    P_c[pcIR.iICIT] = 1 + H_c/Kh[pcIR.iICIT] + Mg_c/Km[pcIR.iICIT]
    P_x[pcIR.iSCOA] = 1 + H_x/Kh[pcIR.iSCOA]
    P_i[pcIR.iSCOA] = 1 + H_i/Kh[pcIR.iSCOA]
    P_c[pcIR.iSCOA] = 1 + H_c/Kh[pcIR.iSCOA]
    P_x[pcIR.iSUC] = 1 + H_x/Kh[pcIR.iSUC] + Mg_x/Km[pcIR.iSUC] + K_x/Kk[pcIR.iSUC]
    P_i[pcIR.iSUC] = 1 + H_i/Kh[pcIR.iSUC] + Mg_i/Km[pcIR.iSUC] + K_i/Kk[pcIR.iSUC]
    P_c[pcIR.iSUC] = 1 + H_c/Kh[pcIR.iSUC] + Mg_c/Km[pcIR.iSUC] + K_c/Kk[pcIR.iSUC]
    P_x[pcIR.iFUM] = 1 + H_x/Kh[pcIR.iFUM]
    P_i[pcIR.iFUM] = 1 + H_i/Kh[pcIR.iFUM]
    P_c[pcIR.iFUM] = 1 + H_c/Kh[pcIR.iFUM]
    P_x[pcIR.iMAL] = 1 + H_x/Kh[pcIR.iMAL] + Mg_x/Km[pcIR.iMAL] + K_x/Kk[pcIR.iMAL]
    P_i[pcIR.iMAL] = 1 + H_i/Kh[pcIR.iMAL] + Mg_i/Km[pcIR.iMAL] + K_i/Kk[pcIR.iMAL]
    P_c[pcIR.iMAL] = 1 + H_c/Kh[pcIR.iMAL] + Mg_c/Km[pcIR.iMAL] + K_c/Kk[pcIR.iMAL]
    P_x[pcIR.iCO2tot] = 1 + H_x/Kh[pcIR.iCO2tot]
    P_i[pcIR.iCO2tot] = 1 + H_i/Kh[pcIR.iCO2tot]
    P_c[pcIR.iCO2tot] = 1 + H_c/Kh[pcIR.iCO2tot]
    P_x[pcIR.iPYR] = 1 + Mg_x/Km[pcIR.iPYR]
    P_i[pcIR.iPYR] = 1 + Mg_i/Km[pcIR.iPYR]
    P_c[pcIR.iPYR] = 1 + Mg_c/Km[pcIR.iPYR]
    P_x[pcIR.iG6P] = 1 + H_x/Kh[pcIR.iG6P]
    P_i[pcIR.iG6P] = 1 + H_i/Kh[pcIR.iG6P]
    P_c[pcIR.iG6P] = 1 + H_c/Kh[pcIR.iG6P]
    P_x[pcIR.iGLU] = 1 + H_x/Kh[pcIR.iGLU] + Mg_x/Km[pcIR.iGLU]
    P_i[pcIR.iGLU] = 1 + H_i/Kh[pcIR.iGLU] + Mg_i/Km[pcIR.iGLU]
    P_c[pcIR.iGLU] = 1 + H_c/Kh[pcIR.iGLU] + Mg_c/Km[pcIR.iGLU]
    P_x[pcIR.iASP] = 1 + H_x/Kh[pcIR.iASP] + Mg_x/Km[pcIR.iASP]
    P_i[pcIR.iASP] = 1 + H_i/Kh[pcIR.iASP] + Mg_i/Km[pcIR.iASP]
    P_c[pcIR.iASP] = 1 + H_c/Kh[pcIR.iASP] + Mg_c/Km[pcIR.iASP]

    ## I.Flux expressions in the TCA cycle
    ## -------------------------------
    # 1. Pyruvate dehydrogenase
    # PYR + COASH + NAD(+H2O) = CO2tot + SCOA + NADH
    # A - PYR; B - COASH; C - NAD; P - CO2tot; Q - SCOA; R - NADH;

    # load concentrations of reactants and products
    A = PYR_x
    B = COASH_x
    C = NAD_x
    P = CO2tot
    Q = ACCOA_x
    R = NADH_x

    # Letters follow the short forms used in the Supplementary Documentation to
    # Wu et al.'s original paper, which is also an excellent source for reading the
    # system of DEs laid out in a more documented form.

    # dG and Keq vlaues
    dGr_pdho = dGf1[pcIR.iCO2tot] + dGf1[pcIR.iACCOA] + dGf1[pcIR.iNADH] \
               - dGf1[pcIR.iPYR] - dGf1[pcIR.iCOASH] - \
               dGf1[pcIR.iNAD] - dGf1[pcIR.iH2O]
    Keq_pdho = np.exp(-dGr_pdho / pcPC.RT)
    Keq_pdh = Keq_pdho * (1 / H_x) * (P_x[pcIR.iCO2tot] *
                        P_x[pcIR.iACCOA] * P_x[pcIR.iNADH])/ (P_x[pcIR.iPYR] *
                        P_x[pcIR.iCOASH] * P_x[pcIR.iNAD])

    # Km and Ki values (Molar)
    KmA = 38.3e-6
    KmB = 9.9e-6
    KmC = 60.7e-6
    KiACCOA = 40.2e-6
    KiNADH = 40.0e-6

    # Inhibition constants
    ai1 = 1 + ACCOA_x / KiACCOA
    ai2 = 1 + NADH_x / (KiNADH * DCA)

    # Vm values
    Vmf = Vmax[pcITCA.ipdh]

    # total reaction flux
    if (A > MinCon) and (B > MinCon) and (C > MinCon):
        J_pdh = Vmf*(A*B*C-P*Q*R/Keq_pdh)/ (KmC*ai2*A*B + KmB*ai1*A*C + KmA*B*C + A*B*C)
    else:
        J_pdh = 0

    # -------------------------------
    # 2. Citrate synthetase
    # OAA + ACCOA(+H2O) = COASH + CIT
    # A - OAA; B - ACCOA; P - COASH; Q - CIT;

    A = OAA_x
    B = ACCOA_x
    P = COASH_x
    Q = CIT_x

    # dG and Keq values
    dGr_citso = dGf1[pcIR.iCOASH] + dGf1[pcIR.iCIT] - dGf1[pcIR.iACCOA] - \
                dGf1[pcIR.iOAA] - dGf1[pcIR.iH2O]
    Keq_citso = np.exp(-dGr_citso / pcPC.RT)
    Keq_cits = Keq_citso * (1 / H_x ** 2) * (P_x[pcIR.iCOASH] * P_x[pcIR.iCIT])\
               / (P_x[pcIR.iACCOA] * P_x[pcIR.iOAA])

    # Km and Ki values (Molar)
    KmA = 4e-6
    KmB = 14e-6
    Kia = 3.33e-6
    KiCIT = 1600e-6
    KiATP = 900e-6
    KiADP = 1800e-6
    KiAMP = 6000e-6
    KiCOASH = 67e-6
    KiSCOA = 140e-6

    # inhibition coefficients
    uCIT_x = CIT_x * (1 + H_x / Kh[pcIR.iCIT]) / P_x[pcIR.iCIT] # unchelated
    uATP_x = ATP_x * (1 + H_x / Kh[pcIR.iATP]) / P_x[pcIR.iATP] # unchelated
    uADP_x = ADP_x * (1 + H_x / Kh[pcIR.iADP]) / P_x[pcIR.iADP] # unchelated
    uAMP_x = AMP_x * (1 + H_x / Kh[pcIR.iAMP]) / P_x[pcIR.iAMP] # unchelated
    ai1 = 1 + uCIT_x / KiCIT
    ai2 = 1 + uATP_x / KiATP + uADP_x / KiADP + uAMP_x / KiAMP\
          + COASH_x / KiCOASH + SCOA_x / KiSCOA

    # Vm values
    Vmf = Vmax[pcITCA.icits]

    # Forward reaction flux
    J_cits_f = Vmf * A * B / (Kia * KmB * ai1 + KmA * ai1 * B + KmB * ai2 * A + A * B)

    # overall reaction flux
    J_cits = J_cits_f - Vmf*(P*Q/Keq_cits) / (Kia*KmB*ai1 + KmA*ai1*B + KmB*ai2*A + A*B)

    # -------------------------------
    # 3. Aconitase
    # CIT = ICIT
    # A - CIT; P - ICIT;

    # load concentrations of reactants and products
    A = CIT_x
    P = ICIT_x

    # dG and Keq values
    dGr_acono = dGf1[pcIR.iICIT] - dGf1[pcIR.iCIT]
    Keq_acono = np.exp(-dGr_acono / pcPC.RT)
    Keq_acon = Keq_acono * P_x[pcIR.iICIT] / P_x[pcIR.iCIT]

    # Km and Ki values (Molar)
    KmA = 1161e-6
    KmP = 434e-6

    # Vm values
    Vmf = Vmax[pcITCA.iacon]
    Vmr = Vmf * (KmP / KmA / Keq_acon)

    # forward reaction flux
    J_acon_f = Vmf*Vmr*A /(KmA*Vmr+Vmr*A+Vmf/Keq_acon*P)

    #  total reaction flux
    J_acon = J_acon_f - Vmf*Vmr*(P/Keq_acon)/(KmA*Vmr+Vmr*A+Vmf/Keq_acon*P)

    # -------------------------------
    # 4. Isocitrate dehydrogenase
    # NAD + ICIT (+ H2O) =  AKG + NADH + CO2tot
    # A - NAD; B - ICIT; P - AKG; Q - NADH; R - CO2tot;

    # load concentrations of reactants and products
    A = NAD_x
    B = ICIT_x
    P = AKG_x
    Q = NADH_x
    R = CO2tot

    # dG and Keq values
    dGr_isodo = dGf1[pcIR.iAKG] + dGf1[pcIR.iNADH] + \
                dGf1[pcIR.iCO2tot] - dGf1[pcIR.iICIT] \
                - dGf1[pcIR.iNAD] - dGf1[pcIR.iH2O]
    Keq_isodo = np.exp(-dGr_isodo / pcPC.RT)
    Keq_isod = Keq_isodo * (1 / H_x ** 2) * (P_x[pcIR.iAKG] *
               P_x[pcIR.iNADH] * P_x[pcIR.iCO2tot])\
               / (P_x[pcIR.iICIT] * P_x[pcIR.iNAD])

    # Km and Ki values (Molar)
    KmA = 74e-6
    KmB = 183e-6
    nH = 3.0
    Kib = 23.8e-6
    Kiq = 29e-6
    KiATP = 91e-6
    KaADP = 50e-6

    # inhibition coefficients
    fATP_x = ATP_x * (1 + H_x / Kh[pcIR.iATP]) / P_x[pcIR.iATP]
    fADP_x = ADP_x * (1 + H_x / Kh[pcIR.iADP]) / P_x[pcIR.iADP]
    ai = 1 + KaADP / fADP_x * (1 + fATP_x / KiATP)

    # Vm values
    Vmf = Vmax[pcITCA.iisod]

    # total reaction flux
    if (A > MinCon) and (B > MinCon):
        J_isod = Vmf / (1 + (KmB / B) ** nH * ai + KmA / A *
                        (1 + (Kib / B) ** nH * ai + Q * ai / Kiq)) * \
                        (1 - 1 / Keq_isod * P * Q * R / A / B)
    else:
        J_isod = 0

    # -------------------------------
    # 5. alpha-Ketoglutarate dehydrogenase
    # AKG + COASH + NAD (+ H2O) = CO2tot + SCOA + NADH
    # A - AKG; B - COASH; C - NAD; P - CO2tot; Q - SCOA; R - NADH;

    # load concentrations of reactants and products
    A = AKG_x
    B = COASH_x
    C = NAD_x
    P = CO2tot
    Q = SCOA_x
    R = NADH_x

    # dG and Keq values
    dGr_akgdo = dGf1[pcIR.iCO2tot] + dGf1[pcIR.iSCOA] + \
                dGf1[pcIR.iNADH] - dGf1[pcIR.iAKG] \
                - dGf1[pcIR.iCOASH] - dGf1[pcIR.iNAD] - dGf1[pcIR.iH2O]
    Keq_akgdo = np.exp(-dGr_akgdo / pcPC.RT)
    Keq_akgd = Keq_akgdo * (1 / H_x) * (P_x[pcIR.iCO2tot] *
               P_x[pcIR.iSCOA] * P_x[pcIR.iNADH])\
               / (P_x[pcIR.iAKG] * P_x[pcIR.iCOASH] * P_x[pcIR.iNAD])

    # Km and Ki values (Molar)
    KmA = 80e-6
    KmB = 55e-6
    KmC = 21e-6
    Kiq = 6.9e-6
    Kir1 = 4.5e-6
    Kir2 = 12.7e-6
    KiATP = 50e-6
    KaADP = 100e-6

    # inhibition coefficients
    fATP_x = ATP_x * (1+H_x/Kh[pcIR.iATP])/P_x[pcIR.iATP]
    fADP_x = ADP_x * (1+H_x/Kh[pcIR.iADP])/P_x[pcIR.iADP]
    ai = 1 + KaADP/fADP_x*(1+fATP_x/KiATP)

    # Vm values
    Vmf = Vmax[pcITCA.iakgd]

    # for testing ## I don't know what they're "testing"
    # Kir1 = param(23)
    # Kir2 = 1e3

    # total reaction flux
    if (A > MinCon) and (B > MinCon) and (C > MinCon):
        J_akgd = Vmf/(1+KmA/A*ai+KmB/B*(1+Q/Kiq)+KmC/C*(1+R/Kir1))/\
                 (1+R/Kir2)*(1-1/Keq_akgd*P*Q*R/A/B/C)
    else:
        J_akgd = 0

    # -------------------------------
    # 6. Succinyl-CoA synthetase
    # GDP + SCOA + PI = COASH + SUC + GTP

    # load concentrations of reactants and products

    A = GDP_x
    B = SCOA_x
    C = Pi_x
    P = COASH_x
    Q = SUC_x
    R = GTP_x

    # dG and Keq values
    dGr_scoaso = dGf1[pcIR.iCOASH] + dGf1[pcIR.iSUC] + dGf1[pcIR.iGTP] \
                 - dGf1[pcIR.iGDP] - dGf1[pcIR.iSCOA] - dGf1[pcIR.iPi]
    Keq_scoaso = np.exp(-dGr_scoaso / pcPC.RT)
    Keq_scoas = Keq_scoaso * (1 / H_x) * (P_x[pcIR.iCOASH] * P_x[pcIR.iSUC]
                * P_x[pcIR.iGTP]) / (P_x[pcIR.iSCOA] * P_x[pcIR.iPi]
                * P_x[pcIR.iGDP])

    # Km and Ki values (Molar)
    KmA = 16e-6
    KmB = 55e-6
    KmC = 660e-6
    KmP = 20e-6
    KmQ = 880e-6
    KmR = 11.1e-6
    Kia = 5.5e-6
    Kib = 100e-6
    Kic = 2000e-6
    Kip = 20e-6
    Kiq = 3000e-6
    Kir = 11.1e-6

    # Vm values
    Vmf = Vmax[pcITCA.iscoas]
    Vmr = Vmf / Keq_scoas * KmP * Kiq * Kir / (Kia * Kib * KmC)

    # total reaction flux
    J_scoas = (Vmf * Vmr * A * B * C - Vmf * Vmr * (P * Q * R / Keq_scoas))/\
              (Vmr * Kia * Kib * KmC + Vmr * Kib * KmC * A + Vmr * Kia * KmB * C
       + Vmr * KmC * A * B + Vmr * KmB * A * C + Vmr * KmA * B * C + Vmr * A * B * C
       + Vmf * Kir * KmQ * P / Keq_scoas + Vmf * Kiq * KmP * R / Keq_scoas + Vmf * KmR
       * P * Q / Keq_scoas + Vmf * KmQ * P * R / Keq_scoas
       + Vmf * KmP * Q * R / Keq_scoas + Vmf * P * Q * R / Keq_scoas + Vmf * KmQ *
       Kir * A * P / Kia / Keq_scoas + Vmr * Kia * KmB * C * R / Kir
       + Vmf * KmQ * Kir * A * B * P / Kia / Kib / Keq_scoas + Vmr * KmA * B * C *
       R / Kir + Vmf * KmR * A * P * Q / Kia / Keq_scoas
       + Vmr * Kia * KmB * C * Q * R / Kiq / Kir + Vmf * Kir * KmQ * A * B * C * P
       / Kia / Kib / Kic / Keq_scoas + Vmf * Kip * KmR * A * B * C * Q / Kia / Kib /
       Kic / Keq_scoas + Vmf * KmR * A * B * P * Q / Kia / Kib / Keq_scoas + Vmr * KmA *
       B * C * Q * R / Kiq / Kir + Vmr * KmA * Kic * B * P * Q * R / Kip / Kiq / Kir
       + Vmr * Kia * KmB * C * P * Q * R / Kip / Kiq / Kir + Vmf * KmR * A * B * C *
       P * Q / Kia / Kib / Kic / Keq_scoas + Vmr * KmA * B * C * P * Q * R / Kip /
       Kiq / Kir)

    # -------------------------------
    # 7. Succinate dehydrogenase
    # SUC + COQ = QH2 + FUM

    # load concentrations of reactants and products
    A = GDP_x
    B = SCOA_x
    C = Pi_x
    P = COASH_x
    Q = SUC_x
    R = GTP_x

    # dG and Keq values
    dGr_scoaso = dGf1[pcIR.iCOASH] + dGf1[pcIR.iSUC] + \
                 dGf1[pcIR.iGTP] - dGf1[pcIR.iGDP] \
                 - dGf1[pcIR.iSCOA] - dGf1[pcIR.iPi]
    Keq_scoaso = np.exp(-dGr_scoaso / pcPC.RT)
    Keq_scoas = Keq_scoaso * (1 / H_x) * (P_x[pcIR.iCOASH] *
                P_x[pcIR.iSUC] * P_x[pcIR.iGTP]) \
                / (P_x[pcIR.iSCOA] * P_x[pcIR.iPi] * P_x[pcIR.iGDP])

    # Km and Ki values (Molar)
    KmA = 16e-6
    KmB = 55e-6
    KmC = 660e-6
    KmP = 20e-6
    KmQ = 880e-6
    KmR = 11.1e-6
    Kia = 5.5e-6
    Kib = 100e-6
    Kic = 2000e-6
    Kip = 20e-6
    Kiq = 3000e-6
    Kir = 11.1e-6

    #Vm values
    Vmf = Vmax[pcITCA.iscoas]
    Vmr = Vmf / Keq_scoas * KmP * Kiq * Kir / (Kia * Kib * KmC)

    # total reaction flux ## I take no responsibility for the unreadability
    J_scoas = (Vmf * Vmr * A * B * C - Vmf * Vmr * (P * Q * R / Keq_scoas)) / \
              (Vmr * Kia * Kib * KmC + Vmr * Kib * KmC * A + Vmr * Kia * KmB * C
       + Vmr * KmC * A * B + Vmr * KmB * A * C + Vmr * KmA * B * C + Vmr * A * B * C
       + Vmf * Kir * KmQ * P / Keq_scoas + Vmf * Kiq * KmP * R / Keq_scoas +
       Vmf * KmR * P * Q / Keq_scoas + Vmf * KmQ * P * R / Keq_scoas
       + Vmf * KmP * Q * R / Keq_scoas + Vmf * P * Q * R / Keq_scoas + Vmf * KmQ *
       Kir * A * P / Kia / Keq_scoas + Vmr * Kia * KmB * C * R / Kir
       + Vmf * KmQ * Kir * A * B * P / Kia / Kib / Keq_scoas + Vmr * KmA * B * C *
       R / Kir + Vmf * KmR * A * P * Q / Kia / Keq_scoas + Vmr * Kia * KmB * C * Q
       * R / Kiq / Kir + Vmf * Kir * KmQ * A * B * C * P / Kia / Kib / Kic / Keq_scoas
       + Vmf * Kip * KmR * A * B * C * Q / Kia / Kib / Kic / Keq_scoas
       + Vmf * KmR * A * B * P * Q / Kia / Kib / Keq_scoas + Vmr * KmA * B * C *
       Q * R / Kiq / Kir + Vmr * KmA * Kic * B * P * Q * R / Kip / Kiq / Kir
       + Vmr * Kia * KmB * C * P * Q * R / Kip / Kiq / Kir + Vmf * KmR * A * B *
       C * P * Q / Kia / Kib / Kic / Keq_scoas + Vmr * KmA * B * C * P * Q * R / Kip
       / Kiq / Kir)

    # -------------------------------
    # 7. Succinate dehydrogenase
    # SUC + COQ = QH2 + FUM

    # load concentrations of reactants and products
    A = SUC_x
    B = COQ_x
    P = QH2_x
    Q = FUM_x

    # dG and Keq values
    dGr_sdho = dGf1[pcIR.iQH2] + dGf1[pcIR.iFUM] - dGf1[pcIR.iSUC] - dGf1[pcIR.iCOQ]
    Keq_sdho = np.exp(-dGr_sdho/pcPC.RT)
    Keq_sdh = Keq_sdho * (P_x[pcIR.iFUM] * P_x[pcIR.iQH2]) / (P_x[pcIR.iSUC] *
                          P_x[pcIR.iCOQ])

    # Km and Ki values(Molar)
    KmA = 467e-6
    KmB = 480e-6
    KmP = 2.45e-6
    KmQ = 1200e-6
    Kia = 120e-6
    Kiq = 1275e-6
    KiOAA = 1.5e-6
    # % % from Gopher and Gutman
    # % KaSUC = 800e-6;
    # % KaFUM = 6400e-6;
    # % from Kohn et al.
    KaSUC = 450e-6
    KaFUM = 375e-6

    # inhibition coefficients
    ai = (1+OAA_x/KiOAA+SUC_x/KaSUC+FUM_x/KaFUM)/(1+SUC_x/KaSUC+FUM_x/KaFUM)

    # Vm values
    Vmf = Vmax[pcITCA.isdh]
    Vmr = Vmf / Keq_sdh * (KmP * Kiq / Kia / KmB)

    # total reaction flux
    J_sdh = (Vmf*Vmr*A*B - Vmf*Vmr*(P*Q/Keq_sdh)) / (Vmr*Kia*KmB*ai+Vmr*KmB*A
                +Vmr*KmA*ai*B+Vmf*KmQ*ai/Keq_sdh*P+Vmf*KmP/Keq_sdh*Q
                +Vmr*A*B+Vmf*KmQ/Kia/Keq_sdh*A*P+Vmr*KmA/Kiq*B*Q
                +Vmf/Keq_sdh*P*Q)

    # -------------------------------
    # 8. Fumarase
    # FUM (+ H2O) = MAL

    # load concentrations of reactants and products
    A = FUM_x
    P = MAL_x

    # dG and Keq values
    dGr_fumo = dGf1[pcIR.iMAL] - dGf1[pcIR.iFUM] - dGf1[pcIR.iH2O]
    Keq_fumo = np.exp(-dGr_fumo/pcPC.RT)
    Keq_fum = Keq_fumo * P_x[pcIR.iMAL] / P_x[pcIR.iFUM]

    # Km and Ki values (Molar)
    KmA = 44.7e-6
    KmP = 197.7e-6
    KiCIT = 3500e-6
    KiATP = 40e-6
    KiADP = 400e-6
    KiGTP = 80e-6
    KiGDP = 330e-6

    # inhibition coefficients
    fATP_x = ATP_x * (1 + H_x / Kh[pcIR.iATP]) / P_x[pcIR.iATP]
    fADP_x = ADP_x * (1 + H_x / Kh[pcIR.iADP]) / P_x[pcIR.iADP]
    fGTP_x = GTP_x * (1 + H_x / Kh[pcIR.iGTP]) / P_x[pcIR.iGTP]
    fGDP_x = GDP_x * (1 + H_x / Kh[pcIR.iGDP]) / P_x[pcIR.iGDP]
    ai = 1 + CIT_x / KiCIT + fATP_x / KiATP + fADP_x / \
         KiADP + fGTP_x / KiGTP + fGDP_x / KiGDP

    # Vm values
    Vmf = Vmax[pcITCA.ifum]
    Vmr = Vmf/Keq_fum*(KmP/KmA)

    # total reaction flux
    J_fum =  (Vmf*Vmr*A - Vmf*Vmr*(P/Keq_fum))/(KmA*Vmr*ai+Vmr*A+Vmf/Keq_fum*P)

    # -------------------------------
    # 9. Malate dehydrogenase
    # NAD + MAL = OAA + NADH(+ H ^ +)

    # load concentrations of reactants and products
    A = NAD_x
    B = MAL_x
    P = OAA_x
    Q = NADH_x

    # dG and Keq values
    dGr_mdho = dGf1[pcIR.iOAA] + dGf1[pcIR.iNADH] - dGf1[pcIR.iNAD] - dGf1[pcIR.iMAL]
    Keq_mdho = np.exp(-dGr_mdho/pcPC.RT)
    Keq_mdh = Keq_mdho * 1 / H_x * (P_x[pcIR.iOAA] * P_x[pcIR.iNADH]) / \
              (P_x[pcIR.iMAL] * P_x[pcIR.iNAD])

    # Km and Ki values (Molar)
    KmA = 90.55e-6
    KmB = 250e-6
    KmP = 6.128e-6
    KmQ = 2.58e-6
    Kia = 279e-6
    Kib = 360e-6
    Kip = 5.5e-6
    Kiq = 3.18e-6
    # from Kohn et al.
    # KiATP = 709.3e-6;
    # KiADP = 383.2e-6;
    # KiAMP = 793.0e-6;
    # from Oza and Shore

    # inhibition coefficients
    fATP_x = ATP_x * (1+H_x/Kh[pcIR.iATP])/P_x[pcIR.iATP]
    fADP_x = ADP_x * (1+H_x/Kh[pcIR.iADP])/P_x[pcIR.iADP]
    fAMP_x = AMP_x * (1+H_x/Kh[pcIR.iAMP])/P_x[pcIR.iAMP]
    ai = 1+fATP_x/KiATP+fADP_x/KiADP+fAMP_x/KiAMP

    # Vm values
    Vmf = Vmax[pcITCA.imdh]
    Vmr = (Vmf/Keq_mdh*(Kiq*KmP/Kia/KmB))

    # total reaction flux
    J_mdh = (Vmf*Vmr*A*B - Vmf*Vmr*(P*Q/Keq_mdh)) / (Vmr*Kia*KmB*ai+Vmr*KmB*A
                +Vmr*KmA*ai*B+Vmf*KmQ*ai/Keq_mdh*P+Vmf*KmP/Keq_mdh*Q
                +Vmr*A*B+Vmf*KmQ/Kia/Keq_mdh*A*P+Vmf/Keq_mdh*P*Q
                +Vmr*KmA/Kiq*B*Q+Vmr/Kip*A*B*P+Vmf/Kib/Keq_mdh*B*P*Q)

    #  -------------------------------
    # 10. Nucleoside diphosphokinase
    # GTP + ADP = GDP + ATP

    # load concentrations of reactants and products
    A = GTP_x
    B = ADP_x
    P = GDP_x
    Q = ATP_x

    #  dG and Keq values
    dGr_ndko = dGf1[pcIR.iGDP] + dGf1[pcIR.iATP] - dGf1[pcIR.iGTP] - dGf1[pcIR.iADP]
    Keq_ndko = np.exp(-dGr_ndko/pcPC.RT)
    Keq_ndk = Keq_ndko * (P_x[pcIR.iGDP] * P_x[pcIR.iATP]) / \
              (P_x[pcIR.iGTP] * P_x[pcIR.iADP])

    # Km and Ki values (Molar)
    KmA = 111e-6
    KmB = 100e-6
    KmP = 260e-6
    KmQ = 278e-6
    Kia = 170e-6
    Kib = 143.6e-6
    Kip = 146.6e-6
    Kiq = 156.5e-6
    KiAMP = 650e-6

    # inhibition coefficients
    fAMP_x = AMP_x * (1+H_x/Kh[pcIR.iAMP])/P_x[pcIR.iAMP]
    ai = 1 + fAMP_x/KiAMP

    # Vm values
    Vmf = Vmax[pcITCA.indk]
    Vmr = Vmf / Keq_ndk * (KmQ * Kip / Kia / KmB)

    # forward reaction flux
    if (A > MinCon) and (B > MinCon):
        J_ndk_f = Vmf*Vmr*A*B /ai / (Vmr*KmB*A+Vmr*KmA*B+Vmf*KmQ/Keq_ndk
                    *P+Vmf*KmP/Keq_ndk*Q+Vmr*A*B+Vmf*KmQ/Kia
                    /Keq_ndk*A*P+Vmf/Keq_ndk*P*Q+Vmr*KmA/Kiq*B*Q)
    else:
        J_ndk_f = 0

    # total reaction flux
    if (P > MinCon) and (Q > MinCon):
        J_ndk = J_ndk_f - Vmf*Vmr*(P*Q/Keq_ndk)/ai / (Vmr*KmB*A+Vmr*KmA*B+Vmf*
                KmQ/Keq_ndk*P+Vmf*KmP/Keq_ndk*Q+Vmr*A*B+Vmf*KmQ/Kia/Keq_ndk*A*P
                +Vmf/Keq_ndk*P*Q+Vmr*KmA/Kiq*B*Q)
    else:
        J_ndk = J_ndk_f

    #  -------------------------------
    #  11. Glutamate oxaloacetate transaminase (aspartate transaminase)
    #  ASP + AKG = OAA + GLU

    # load concentrations of reactants and products
    A = ASP_x
    B = AKG_x
    P = OAA_x
    Q = GLU_x

    # dG and Keq values
    dGr_goto = dGf1[pcIR.iOAA] + dGf1[pcIR.iGLU] - dGf1[pcIR.iASP] - dGf1[pcIR.iAKG]
    Keq_goto = np.exp(-dGr_goto / pcPC.RT)
    Keq_got = Keq_goto * (P_x[pcIR.iOAA] * P_x[pcIR.iGLU]) / (P_x[pcIR.iASP] *
                                                              P_x[pcIR.iAKG])

    #  Km and Ki values (Molar)
    KmA = 3900e-6
    KmB = 430e-6
    KmP = 88e-6
    KmQ = 8900e-6
    Kia = 3480e-6
    Kib = 710e-6
    Kip = 50e-6
    Kiq = 8400e-6
    KiAKG = 16.6e-3
    ai = 1 + AKG_x / KiAKG

    # Vm values
    Vmf = Vmax[pcITCA.igot]
    Vmr = Vmf/Keq_got*(KmQ*Kip/Kia/KmB)

    # forward reaction flux
    if (A > MinCon) and (B > MinCon):
        J_got_f = Vmf*Vmr*A*B / (Vmr*KmB*A+Vmr*KmA*ai*B
                +Vmf*KmQ/Keq_got*ai*P+Vmf*KmP/Keq_got*Q+Vmr*A*B
                +Vmf*KmQ/Kia/Keq_got*A*P+Vmf/Keq_got*P*Q
                +Vmr*KmA/Kiq*B*Q)
    else:
        J_got_f = 0

    # total reaction flux
    if (P > MinCon) and (Q > MinCon):
        J_got = J_got_f - Vmf*Vmr*(P*Q/Keq_got) / (Vmr*KmB*A+Vmr*KmA*ai*B
                +Vmf*KmQ/Keq_got*ai*P+Vmf*KmP/Keq_got*Q+Vmr*A*B
                +Vmf*KmQ/Kia/Keq_got*A*P+Vmf/Keq_got*P*Q
                +Vmr*KmA/Kiq*B*Q)
    else:
        J_got = J_got_f

    # -------------------------------
    # 12. Anti- and co-transporter fluxes of substrates involved in the TCA cycle
    # Fluxes are defined to be positive when the first reactant(s) move from IM
    # into mito. matrix

    # --------------------------------
    # (1) PYR^{-}-H^{+} co-transporter

    PYR_i1 = PYR_i * 1 / P_i[pcIR.iPYR]
    PYR_x1 = PYR_x * 1 / P_x[pcIR.iPYR]
    J_PYR_H = x_PYR_H * (PYR_i1 * H_i - PYR_x1 * H_x)

    # --------------------------------
    # (2) GLU^{-}-H^{+} co-transporter

    GLU_i1 = GLU_i * 1 / P_i[pcIR.iGLU]
    GLU_x1 = GLU_x * 1 / P_x[pcIR.iGLU]
    J_GLU_H = x_GLU_H * (GLU_i1 * H_i - GLU_x1 * H_x)

    # --------------------------------
    # (3) CIT^{2-}/MAL^{2-} anti-transporter

    CIT_i1 = CIT_i * (H_i / Kh[pcIR.iCIT]) / P_i[pcIR.iCIT]
    CIT_x1 = CIT_x * (H_x / Kh[pcIR.iCIT]) / P_x[pcIR.iCIT]
    MAL_i1 = MAL_i * 1 / P_i[pcIR.iMAL]
    MAL_x1 = MAL_x * 1 / P_x[pcIR.iMAL]
    J_CIT_MAL = x_CIT_MAL * (CIT_i1 * MAL_x1 - CIT_x1 * MAL_i1)

    # --------------------------------
    # (4) AKG^{2-}/MAL^{2-} anti-transporter
    KmMALi = 1.4e-3
    KmMALx = 0.7e-3
    KmAKGi = 0.3e-3
    KmAKGx = 0.17e-3

    J_AKG_MAL = x_AKG_MAL/(KmMALx*KmAKGi)*(MAL_x*AKG_i-MAL_i*AKG_x)/\
         (2 + MAL_i/KmMALi + MAL_x/KmMALx + AKG_i/KmAKGi + AKG_x/KmAKGx +
          MAL_i*AKG_x/(KmMALi*KmAKGx) + MAL_x*AKG_i/(KmMALx*KmAKGi))

    # --------------------------------
    # (5) MAL^{2-}/PI^{2-} anti-transporter

    MAL_i1 = MAL_i*1/P_i[pcIR.iMAL]
    MAL_x1 = MAL_x*1/P_x[pcIR.iMAL]
    Pi_i1 = Pi_i*1/P_i[pcIR.iPi]
    Pi_x1 = Pi_x*1/P_x[pcIR.iPi]
    J_MAL_PI = x_MAL_PI * (MAL_i1*Pi_x1 - MAL_x1*Pi_i1)

    # --------------------------------
    # (6) ASP^{-}/H.GLU{0} antitransporter (with one negative charge
    # translocated from IM into mito matrix)
    Kiaspi = 28e-6
    Kiaspx = 2.8e-3
    Kiglui = 180e-6
    Kiglux = 1.6e-3
    pKa_gaa = 6.5
    Kh_gaa = 10.**(-pKa_gaa)
    Keq_gaa = np.exp(-pcPC.F*dPsi/pcPC.RT)*P_x[pcIR.iASP]*\
              P_i[pcIR.iGLU]/P_i[pcIR.iASP]/P_x[pcIR.iGLU]
    m = 1.8
    J_ASP_GLU = x_ASP_GLU/(Keq_gaa*Kiaspi*Kiglux*Kh_gaa)*(Keq_gaa*ASP_i*GLU_x*
            H_x - ASP_x*GLU_i*H_i)/ (2*m + m*ASP_i/Kiaspi + ASP_i*GLU_x*H_x/(Kiaspi
            *Kiglux*Kh_gaa) + m*ASP_x*H_i/(Kiaspx*Kh_gaa) + ASP_x*GLU_i*H_i/(Kiaspx*
            Kiglui*Kh_gaa) + m*ASP_x/Kiaspx + m*ASP_i*H_x/(Kiaspi*Kh_gaa) +
            m*H_x/Kh_gaa + m*GLU_i*H_i/(Kiglui*Kh_gaa) + m*H_i/Kh_gaa + m*GLU_x*
            H_x/(Kiglux*Kh_gaa))

    #  --------------------------------
    #  (7) SUC^{2-}/MAL^{2-} anti-transporter

    MAL_i1 = MAL_i * 1 / P_i[pcIR.iMAL]
    MAL_x1 = MAL_x * 1 / P_x[pcIR.iMAL]
    SUC_i1 = SUC_i * 1 / P_i[pcIR.iSUC]
    SUC_x1 = SUC_x * 1 / P_x[pcIR.iSUC]
    J_SUC_MAL = x_SUC_MAL * (SUC_i1 * MAL_x1 - SUC_x1 * MAL_i1)

    # -------------------------------
    # 13. Passive permeation between cytoplasm/buffer and IM
    # Fluxes are defined to be positive when the reactant moves from
    # cytoplasm/buffer into IM

    J_PYRt = gamma * x_PYRt * (PYR_c - PYR_i)
    J_CITt = gamma * x_CITt * (CIT_c - CIT_i)
    J_MALt = gamma * x_MALt * (MAL_c - MAL_i)
    J_AKGt = gamma * x_AKGt * (AKG_c - AKG_i)
    J_SUCt = gamma * x_SUCt * (SUC_c - SUC_i)
    J_GLUt = gamma * x_GLUt * (GLU_c - GLU_i)
    J_ASPt = gamma * x_ASPt * (ASP_c - ASP_i)
    J_FUMt = gamma * x_FUMt * (FUM_c - FUM_i)
    J_ICITt = gamma * x_ICITt * (ICIT_c - ICIT_i)

    # -------------------------------
    # 14. Glutamate Dehydrogenase
    #
    # Glutamate consumption producing AKG

    J_gdh = 2e-10*GLU_x*(NAD_x/9.35e-3)/(3e-3+GLU_x)
    if NAD_x < 1e-10 or GLU_x < 1e-10:
        J_gdh = 0.
    #J_gdh = 0.0

    #  II. Flux expresssions in the oxidative phosphorylation

    # -------------------------------
    # Substrate/ion transport
    # Transport of ADP, ATP, and Pi across outer membrane:

    J_ADP = gamma * x_A * (ADP_c - ADP_i)
    J_ATP = gamma * x_A * (ATP_c - ATP_i)
    J_AMP = gamma * x_A * (AMP_c - AMP_i)
    J_Pi2 = gamma * x_Pi2 * (Pi_c - Pi_i)

    # -------------------------------
    # 1. Complex I
    # NADH_x + Q_x + 5H+_x <-> NAD+_x + QH2_x + 4H+_i + 4dPsi
    # compute free energy of reaction from free energe of formation
    dGr_C1o = dGf1[pcIR.iNAD] + dGf1[pcIR.iQH2] - dGf1[pcIR.iNADH] - dGf1[pcIR.iCOQ]
    # compute Keq from combined free energy of reactions(including potential
    # change due to charge translocation)
    Keq_C1 = np.exp(-(dGr_C1o + 4 * pcPC.F * dPsi) / pcPC.RT)
    Kapp_C1 = Keq_C1 * H_x**5 / H_i**4
    J_C1 = w[0] * x_C1 * (Kapp_C1 * NADH_x * COQ_x - NAD_x * QH2_x)

    # -------------------------------
    # 2. Complex III
    # QH2_x + 2cytoC(ox)3+_i + 2H+_x <-> Q_x + 2cytoC(red)2+_i + 4H+_i + 2dPsi

    dGr_C3o = dGf1[pcIR.iCOQ] + 2 * dGf1[pcIR.iCred] - dGf1[pcIR.iQH2] - \
              2 * dGf1[pcIR.iCox]
    Keq_C3 = np.exp(-(dGr_C3o + 2 * pcPC.F * dPsi) / pcPC.RT)
    Kapp_C3 = Keq_C3 * H_x ** 2 / H_i ** 4
    QH2_x = max(MinCon, QH2_x)
    COQ_x = max(MinCon, COQ_x)
    J_C3 = w[1] * x_C3 * ((1 + Pi_x / k_Pi1) / (1 + Pi_x / k_Pi2)) * (Kapp_C3 ** 0.5 *
                                    Cox_i * np.sqrt(QH2_x) - Cred_i * np.sqrt(COQ_x))

    # -------------------------------
    # 3. Complex IV
    # 2 cytoC(red) 2 + _i + 0.5 O2_x + 4 H + _x <-> 2 cytoC(ox) 3 + _x +
    # H2O_x + 2 H + _i + 2 dPsi

    dGr_C4o = 2 * dGf1[pcIR.iCox] + dGf1[pcIR.iH2O] - 2 * dGf1[pcIR.iCred]\
              - 0.5 * dGf1[pcIR.iO2]

    # 2 charges from translocation of proton, and the other 2 from cytoC

    Keq_C4 = np.exp(-(dGr_C4o + 4 * pcPC.F * dPsi) / pcPC.RT)
    Kapp_C4 = Keq_C4 * H_x ** 4 / H_i ** 2
    O2 = max(O2, MinCon)
    J_C4 = w[2] * x_C4 * (O2 / (O2 + pcPC.k_O2)) * np.exp(pcPC.F * dPsi /
                pcPC.RT) * (Cred_i / pcPC.Ctot) * (
                Kapp_C4 ** 0.5 * Cred_i * (O2 ** 0.25) - Cox_i)

    #  -------------------------------
    #  4. F1Fo-ATPase
    #  ADP3-_x + HPO42-_x + H+_x + n_A*H+_i <-> ATP4- + H2O + n_A*H+_x

    dGr_F1o = dGf1[pcIR.iATP] + dGf1[pcIR.iH2O] - dGf1[pcIR.iADP] - dGf1[pcIR.iPi]
    Keq_F1 = np.exp(-(dGr_F1o - pcPC.n_A * pcPC.F * dPsi) / pcPC.RT)
    Kapp_F1 = Keq_F1 * H_i ** pcPC.n_A / H_x ** (pcPC.n_A - 1) * \
              P_x[pcIR.iATP] / (P_x[pcIR.iADP] * P_x[pcIR.iPi])
    J_F1 = w[3] * x_F1 * (Kapp_F1 * ADP_x * Pi_x - ATP_x)

    #  -------------------------------
    #  5. ANT
    #  ATP4-_x + ADP3-_i <-> ATP4-_i + ADP3-_x

    ADP_i1 = ADP_i / P_i[pcIR.iADP] # ADP ^ 3 -
    ATP_i1 = ATP_i / P_i[pcIR.iATP] # ATP ^ 4 -
    ADP_x1 = ADP_x / P_x[pcIR.iADP] # ADP ^ 3 -
    ATP_x1 = ATP_x / P_x[pcIR.iATP] # ATP ^ 4 -

    theta = 0.60
    Psi_i = theta * dPsi
    Psi_x = (theta - 1) * dPsi
    if (ADP_i > 1e-9) and (ATP_x > 1e-9):
        J_ANT = x_ANT * (ADP_i1 / (ADP_i1 + ATP_i1 *
                        np.exp(-pcPC.F * Psi_i / pcPC.RT)) - ADP_x1 /
                        (ADP_x1 + ATP_x1 * np.exp(-pcPC.F * Psi_x /
                        pcPC.RT))) * (ADP_i1 / (ADP_i1 + pcPC.k_mADP))
    else:
        J_ANT = 0

    #  -------------------------------
    #  6. H+-Pi2 cotransporter

    H2PIi1 = Pi_i * (H_i / Kh[pcIR.iPi]) / P_i[pcIR.iPi]
    H2PIx1 = Pi_x * (H_x / Kh[pcIR.iPi]) / P_x[pcIR.iPi]
    J_Pi1 = x_Pi1 / k_PiH * (H_i * H2PIi1 - H_x * H2PIx1) / (1 + H2PIi1 / k_PiH) \
            / (1 + H2PIx1 / k_PiH)

    #  -------------------------------
    #  7a. H+ leak

    if abs(dPsi) > 1e-9:
        J_Hle2 = x_Hle * dPsi * (H_i * np.exp(pcPC.F * dPsi / pcPC.RT) - H_x)\
                 / (np.exp(pcPC.F * dPsi / pcPC.RT) - 1)
    else:
        J_Hle2 = x_Hle * pcPC.RT * (H_i - H_x) / pcPC.F

    if abs(dPsi) > 1e-9:
        J_Hle = P_Hle * (H_i * np.exp(pcPC.F * dPsi / pcPC.RT / 2.0) - H_x *
                         np.exp(-pcPC.F * dPsi / pcPC.RT / 2.0))
    else:
        J_Hle = P_Hle * (H_i - H_x)

    # -------------------------------
    # 7b. K+ leak

    if abs(dPsi) > 1e-9:
        J_Kle2 = x_Kle * dPsi * (K_i * np.exp(pcPC.F * dPsi / pcPC.RT) - K_x)\
                 / (np.exp(pcPC.F * dPsi / pcPC.RT) - 1)
    else:
        J_Kle2 = x_Kle * pcPC.RT * (K_i - K_x) / pcPC.F

    if abs(dPsi) > 1e-9:
        J_Kle = P_Kle * (K_i * np.exp(pcPC.F * dPsi / pcPC.RT / 2.0) -
                         K_x * np.exp(-pcPC.F * dPsi / pcPC.RT / 2.0))
    else:
        J_Kle = P_Kle * (K_i - K_x)

    #  -------------------------------
    #  8. K+/H+ anti-porter

    J_KH   = x_KH*(K_i*H_x - K_x*H_i)

    #  -------------------------------
    #  9. Rate of cytosolic ATP consumption by Hexokinase reaction
    #  GLC + ATP = ADP + G6P
    #  GLC^0 + ATP^4- = ADP^3- + G6P^2- + H^+
    #  Keq_HK = 880; % http://web.mit.edu/esgbio/www/eb/chem/solvingequilibria.html

    dGr_HKo = dGf1[pcIR.iG6P] + dGf1[pcIR.iADP] - dGf1[pcIR.iGLC] - dGf1[pcIR.iATP]
    Keq_HK = np.exp(-(dGr_HKo) / pcPC.RT)
    Kapp_HK = Keq_HK / H_c * P_c[pcIR.iG6P] * P_c[pcIR.iADP]\
              / (P_c[pcIR.iGLC] * P_c[pcIR.iATP])

    KmA = 1e-3 # MgATP
    Kia = 1e-3
    KmB = 47e-6 # GLC
    Kib = 47e-6
    KmP = 47e-6 # G6P
    Kip = 47e-6
    KmQ = 1e-3 # MgADP
    Kiq = 1e-3
    KiG6P = 10e-6 # inhibition constant of G6P
    A = ATP_c * (Mg_c / Km[pcIR.iATP]) / P_c[pcIR.iATP]
    B = GLC_c
    P = G6P_c
    Q = ADP_c * (Mg_c / Km[pcIR.iADP]) / P_c[pcIR.iADP]
    Kapp_HK_m = Kapp_HK * (Km[pcIR.iATP] * P_c[pcIR.iATP])\
                / (Km[pcIR.iADP] * P_c[pcIR.iADP])
    J_HK = x_HK / (Kib * KmA) * (A * B - P * Q / Kapp_HK_m)\
           / (1 + A / Kia + B / Kib + A * B / Kib / KmA + P
              / Kip + Q / Kiq + P * Q / Kiq / KmP + P * B / KiG6P / Kib)

    #  ---------------------------------------
    #  10. Creatine kinase reaction
    #  ADP3- + PCr- + H+ = ATP4- + Cr0
    #  set Cr and PCr concentrations according to experiments

    if (ExpType == 2) or (ExpType == 5):
        J_CKe = 0
    else:
        x_CK = 1e7
        K_CK = np.exp(50.78 / pcPC.RT)
        CRtot = 42.7e-3
        PCr_c = x[pcIS.iPCr_c]
        Cr_c = CRtot - PCr_c
        ATP_c1 = ATP_c * 1 / P_c[pcIR.iATP] # Mg2 + unbound species
        ADP_c1 = ADP_c * 1 / P_c[pcIR.iADP] # Mg2 + unbound species
        J_CKe = x_CK * (K_CK * ADP_c1 * PCr_c * H_c - ATP_c1 * Cr_c)

    # ------------------------------------------
    # 11. Adenylate kinase reaction
    # 2ADP3- = ATP4- + AMP2-

    dGr_AKo = dGf1[pcIR.iATP] + dGf1[pcIR.iAMP] - 2 * dGf1[pcIR.iADP]
    Keq_AK = np.exp(-dGr_AKo / pcPC.RT)
    Kapp_AKi = Keq_AK * P_i[pcIR.iATP] * P_i[pcIR.iAMP] / \
               (P_i[pcIR.iADP] * P_i[pcIR.iADP])
    Kapp_AKc = Keq_AK * P_c[pcIR.iATP] * P_c[pcIR.iAMP] / \
               (P_c[pcIR.iADP] * P_c[pcIR.iADP])

    if (ExpType == 2) or (ExpType == 5):
        J_AKi = 0
        J_AKe = 0
    else:
        x_AK = 1e7
        J_AKi = x_AK * (Kapp_AKi * ADP_i * ADP_i - AMP_i * ATP_i)
        J_AKe = x_AK * (Kapp_AKc * ADP_c * ADP_c - AMP_c * ATP_c)


    # Outputting fluxes
    J = np.zeros(pc.N_flux)
    J[pcFl.J_C1] = J_C1
    J[pcFl.J_C3] = J_C3
    J[pcFl.J_C4] = J_C4
    J[pcFl.J_F1] = J_F1
    J[pcFl.J_ANT] = J_ANT
    J[pcFl.J_Pi1] = J_Pi1
    J[pcFl.J_Hle] = J_Hle
    J[pcFl.J_KH] = J_KH
    J[pcFl.J_pdh] = J_pdh
    J[pcFl.J_cits] = J_cits
    J[pcFl.J_acon] = J_acon
    J[pcFl.J_isod] = J_isod
    J[pcFl.J_akgd] = J_akgd
    J[pcFl.J_scoas] = J_scoas
    J[pcFl.J_sdh] = J_sdh
    J[pcFl.J_fum] = J_fum
    J[pcFl.J_mdh] = J_mdh
    J[pcFl.J_ndk] = J_ndk
    J[pcFl.J_got] = J_got
    J[pcFl.J_Pi2] = J_Pi2
    J[pcFl.J_ATP] = J_ATP
    J[pcFl.J_ADP] = J_ADP
    J[pcFl.J_AMP] = 0. # J_AMP
    J[pcFl.J_PYR_H] = J_PYR_H
    J[pcFl.J_GLU_H] = J_GLU_H
    J[pcFl.J_CIT_MAL] = J_CIT_MAL
    J[pcFl.J_SUC_MAL] = J_SUC_MAL
    J[pcFl.J_AKG_MAL] = J_AKG_MAL
    J[28] = 0.
    J[29] = 0.
    J[pcFl.J_MAL_PI] = J_MAL_PI
    J[pcFl.J_ASP_GLU] = J_ASP_GLU
    J[pcFl.J_PYRt] = J_PYRt
    J[pcFl.J_CITt] = J_CITt
    J[pcFl.J_ICITt] = J_ICITt
    J[pcFl.J_MALt] = J_MALt
    J[pcFl.J_AKGt] = J_AKGt
    J[pcFl.J_SUCt] = J_SUCt
    J[pcFl.J_FUMt] = J_FUMt
    J[pcFl.J_GLUt] = J_GLUt
    J[pcFl.J_ASPt] = J_ASPt
    J[pcFl.J_HK] = J_HK
    J[pcFl.J_Kle] = J_Kle
    J[pcFl.J_AKi] = J_AKi
    J[pcFl.J_CKe] = J_CKe
    J[pcFl.J_AKe] = J_AKe
    J[pcFl.J_gdh] = J_gdh

    return J