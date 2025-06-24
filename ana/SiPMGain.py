import numpy as np
from scipy.optimize import curve_fit

import warnings

def gauss(V, A, mu, sigma):
    return A * np.exp(-(V-mu)**2/(2*sigma**2))

def ngauss(V, *p):
    assert(len(p) % 3 == 0)
    ngaus = len(p) // 3
    return np.sum([gauss(V, *p[i*3:(i+1)*3]) for i in range(ngaus)], axis=0)

def SiPMGain(sipm_data):
    bins = np.linspace(0, 0.05, 101)
    centers = (bins[:-1] + bins[1:]) / 2

    Ns = []
    for i_sipm in range(sipm_data["hit_amp"].shape[0]):
        N,_ = np.histogram(sipm_data["hit_amp"][i_sipm], bins=bins)
        Ns.append(N)

    NFIT = 5
    gains = []
    for i_sipm in range(len(Ns)):
        # fit parameters
        mu_p0 = np.array([0.005*(i+1) for i in range(NFIT)])
        sig_p0 = np.array([0.0005]*NFIT)
        A_p0 = np.array([np.max(Ns[i_sipm])/1.5**i for i in range(NFIT)])
        
        p0 = sum([[A_p0[i], mu_p0[i], sig_p0[i]] for i in range(NFIT)], [])
        
        mu_lo = mu_p0 - 0.0025
        mu_hi = mu_p0 + 0.0025
        
        sig_lo = np.array([0.00025]*NFIT)
        sig_hi = np.array([0.00100]*NFIT)
        
        A_lo = A_p0/2
        A_hi = A_p0*2
        
        lo = sum([[A_lo[i], mu_lo[i], sig_lo[i]] for i in range(NFIT)], [])
        hi = sum([[A_hi[i], mu_hi[i], sig_hi[i]] for i in range(NFIT)], [])

        # Do the fit
        where_fit = (centers > 0.003) & (centers < 0.027)
        popt, perr = curve_fit(ngauss, centers[where_fit], Ns[i_sipm][where_fit], p0=p0, bounds=(lo, hi))
        
        peaks = popt[1::3]
    
        npe = np.arange(1, len(peaks)+1)
        m, b = np.polyfit(npe, peaks, 1)

        gains.append(m)

    return gains

if __name__ == "__main__":
    import os
    from DataHandling import ReadBinary, GetSBCEvent
    from SiPMPulses import SiPMPulses

    path = "/bluearc/storage/SBC-22-caendata"
    file = "202206201045.bin"
    
    filename = os.path.join(path, file)
    
    data = ReadBinary.ReadBlock(filename)
    sipm_data = SiPMPulses(data)
    print(SiPMGain(sipm_data))
