import numpy as np
from solve_ecc import solve_ecc
from gwb_funcs import F,g0,gm,gp
from newtonian_gwb import HH as HHN

def HH(M, eta, f0, e0, f, z, pmax):
    """  
    Returns the sky-averaged contribution of a single source over its entire evolutionary history.
    Has units of s^3 (volume).
    Averaging this quantity over the parameter space using the SMBHB number distribution gives hc^2(f). 
    """
    Mch = M * eta**(3/5)
    xi = 2*np.pi*Mch*f
    return 2/3/np.pi**2/((1+z)**(1/3))/f**3 * xi**(5/3) * Hsum(M, f0, e0, f, z, pmax)

def compute_fp_ep(M, f0, e0, f, z, p, q):
    fr = f/(1+z)
    
    k = 0
    for i in range(5):
        fp = fr/(p+q*k)
        ep = solve_ecc(f0, e0, fp)[0]
        xi = (2*np.pi*M*fp)
        k  = 3*xi**(2/3)/(1-ep**2)
    if k>0.5:
        print(M, f0, e0, f, z, p, q)
    return fp, ep
    
def y_test(M, fp, ep):
    xi = (2*np.pi*M*fp)
    k  = 3*xi**(2/3)/(1-ep**2)
    y = (M*(1+k)*2*np.pi*fp)**(1/3) / np.sqrt(1-ep**2)
    #if y>0.1:
    #    print(M, fp, ep, y)
    return int(k<0.5)

def Hsum(M, f0, e0, f, z, pmax):
    result = 0
    
    for p in range(1,pmax+1):
        fp,ep = compute_fp_ep(M, f0, e0, f, z, p, -2)
        
        Hterm = gm(p,ep)/F(ep)/p**(2/3)
        result += y_test(M, fp, ep) * Hterm

        fp,ep = compute_fp_ep(M, f0, e0, f, z, p, 0)
        Hterm = g0(p,ep)/F(ep)/p**(2/3)
        result += y_test(M, fp, ep) * Hterm
        
        fp,ep = compute_fp_ep(M, f0, e0, f, z, p, 2)
        Hterm = gp(p,ep)/F(ep)/p**(2/3)
        result += y_test(M, fp, ep) * Hterm
        
    return result
    
if __name__=='__main__':
    import matplotlib.pyplot as plt
    
    Msun = 4.92703806e-6 # s
    Mpc = 1.02927125e14 # s
    
    Ms = np.array([1e7, 1e8, 1e9])*Msun
    eta = 0.25
    f0 = 1e-9 # Hz
    e0s = [0.1, 0.5, 0.9]
    
    fs = 10**np.linspace(-10,-6.5, 20)
    
    z = 0.01
    pmax = 50
    
    iplt = 1
    for M in Ms:
        for e0 in e0s:
            print(M,e0)
            HHs = np.array([HH(M, eta, f0, e0, f, z, pmax) for f in fs])
            HHNs = np.array([HHN(M, eta, f0, e0, f, z, pmax) for f in fs])
            
            plt.subplot(3,3,iplt)
            label="M={:.0e} MSun, e0={:.1f}".format(M/Msun,e0)
            plt.loglog(fs, np.sqrt(HHs/Mpc**3), label=label)
            plt.loglog(fs, np.sqrt(HHNs/Mpc**3), label=label+"(N)")
            plt.grid()
            plt.legend()
            
            if e0==e0s[0]:
                plt.ylabel("$\\sqrt{\\mathcal{H}(f)}$ (Mpc$^{3/2}$)")
            if M==Ms[-1]:
                plt.xlabel("$f$ (Hz)")
            
            iplt += 1
    
    plt.show()
