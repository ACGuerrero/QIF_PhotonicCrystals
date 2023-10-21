import numpy as np
from scipy.constants import c
import matplotlib.pyplot as plt

def wavelength_to_pulsation(wl):
    return (2*np.pi*c)/(wl)

def pulsation_to_wavelength(om):
    return (2*np.pi*c)/(om)

def get_max_bandgap(n1,n2,a,target_puls):
    wm = 2*np.pi*c*(N1+N2)/(a*4*N1*N2)
    l = 2*np.pi*c/wm
    indzero = np.where((l>(target_central_wl-epsilon)) & (l<(target_central_wl+epsilon))==False)
    mask = np.ones(N1.shape)
    mask[indzero] = 0
    Dw =(4/(np.pi))*np.arcsin(np.abs(N1-N2)/(N1+N2))*target_puls
    ind_for_disp = np.where(Dw*mask==np.max(Dw*mask))
    return (pulsation_to_wavelength(np.max(Dw*mask)), n1[ind_for_disp],n2[ind_for_disp])