import numpy as np
from gwb_funcs import F,g
from solve_ecc import solve_ecc

def Hsum(f0, e0, f, z, pmax):
    result = 0
    
    for p in range(1,pmax+1):
        fp = f/(1+z)/p
        ep = solve_ecc(f0, e0, fp)[0]
        Hterm = g(p,ep)/F(ep)/p**(2/3)
        result += Hterm
    
    return result

def HH(M, eta, f0, e0, f, z, pmax):
    """  
    Returns the sky-averaged contribution of a single source over its entire evolutionary history.
    Has units of s^3 (volume).
    Averaging this quantity over the parameter space using the SMBHB number distribution gives hc^2(f). 
    """
    Mch = M * eta**(3/5)
    xi = 2*np.pi*Mch*f
    return 2/3/np.pi**2/((1+z)**(1/3))/f**3 * xi**(5/3) * Hsum(f0, e0, f, z, pmax)
    
    
if __name__=='__main__':
    import matplotlib.pyplot as plt
    
    Msun = 4.92703806e-6 # s
    Mpc = 1.02927125e14 # s
    
    Mchs = np.array([1e7, 1e8, 1e9])*Msun
    f0 = 1e-9 # Hz
    e0s = [0.1, 0.5, 0.9]
    
    fs = 10**np.linspace(-10,-6,50)
    
    z = 0.01
    pmax = 200
    
    iplt = 1
    for Mch in Mchs:
        for e0 in e0s:
            HHs = np.array([HH(Mch, f0, e0, f, z, pmax) for f in fs])
            
            plt.subplot(330+iplt)
            label="Mch={:.0e} MSun, e0={:.1f}".format(Mch/Msun,e0)
            plt.loglog(fs, np.sqrt(HHs/Mpc**3), label=label)
            plt.grid()
            plt.legend()
            
            if e0==e0s[0]:
                plt.ylabel("$\\sqrt{\\mathcal{H}(f)}$ (Mpc$^{3/2}$)")
            if Mch==Mchs[-1]:
                plt.xlabel("$f$ (Hz)")
            
            iplt += 1
    
    plt.show()
    
    
