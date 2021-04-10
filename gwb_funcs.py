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
                     

