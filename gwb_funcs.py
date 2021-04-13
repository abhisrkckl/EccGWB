import numpy as np
from scipy.special import jv as BesselJ

def F(e):
    return (1+73/34*e**2+37/96*e**4) / (1-e**2)**3.5
    
def g(p,e):
    pe = p*e
    Jpm2_pe = BesselJ(p-2, pe)
    Jpm1_pe = BesselJ(p-1, pe)
    Jp_pe   = BesselJ(p,   pe)
    Jpp1_pe = BesselJ(p+1, pe)
    Jpp2_pe = BesselJ(p+2, pe)

    return p**4/32 * (              (Jpm2_pe - 2*e*Jpm1_pe + 2/p*Jp_pe + 2*e*Jpp1_pe - Jpp2_pe)**2
                        + (1-e**2)* (Jpm2_pe               - 2*  Jp_pe               + Jpp2_pe)**2
                        + 4/3/p**2* (                            Jp_pe                        )**2 
                     )
                     
def gp(p,e):
    pe = p*e
    OTS = np.sqrt(1-e**2)
    Jpm2_pe = BesselJ(p-2, pe)
    Jpm1_pe = BesselJ(p-1, pe)
    Jp_pe   = BesselJ(p,   pe)
    Jpp1_pe = BesselJ(p+1, pe)
    Jpp2_pe = BesselJ(p+2, pe)
    
    gp = p**4/64 * (      Jpm2_pe * (1 + OTS)
                        + Jpm1_pe * (-2*e)
                        + Jp_pe   * (2/p - 2*OTS)
                        + Jpp1_pe * (2*e)
                        + Jpp2_pe * (-1 + OTS)  )**2
    
    return gp

def g0(p,e):
    pe = p*e
    #OTS = np.sqrt(1-e**2)
    #Jpm2_pe = BesselJ(p-2, pe)
    #Jpm1_pe = BesselJ(p-1, pe)
    Jp_pe   = BesselJ(p,   pe)
    #Jpp1_pe = BesselJ(p+1, pe)
    #Jpp2_pe = BesselJ(p+2, pe)
   
    g0 = p**2/24 * (      Jp_pe                 )**2
    
    return g0
    
def gm(p,e):
    pe = p*e
    OTS = np.sqrt(1-e**2)
    Jpm2_pe = BesselJ(p-2, pe)
    Jpm1_pe = BesselJ(p-1, pe)
    Jp_pe   = BesselJ(p,   pe)
    Jpp1_pe = BesselJ(p+1, pe)
    Jpp2_pe = BesselJ(p+2, pe)

    gm = p**4/64 * (      Jpm2_pe * (1 - OTS)
                        + Jpm1_pe * (-2*e)
                        + Jp_pe   * (2/p + 2*OTS)
                        + Jpp1_pe * (2*e)
                        + Jpp2_pe * (-1 - OTS)  )**2
   
    return gm
