import numpy as np
from scipy.constants import c
import matplotlib.pyplot as plt

def wavelength_to_pulsation(wl):
    return (2*np.pi*c)/(wl)

def pulsation_to_wavelength(om):
    return (2*np.pi*c)/(om)

def get_max_bandgap(n1_given,n2_given,a,target_wl_user):
    wm = 2*np.pi*c*(n1_given+n2_given)/(a*4*n1_given*n2_given)
    l = 2*np.pi*c/wm
    indzero = np.where((l>(target_wl_user-epsilon)) & (l<(target_wl_user+epsilon))==False)
    mask = np.ones(N1.shape)
    mask[indzero] = 0
    Dw =(4/(np.pi))*np.arcsin(np.abs(n1_given-n2_given)/(n1_given+n2_given))*wavelength_to_pulsation(target_wl_user)
    ind_for_disp = np.where(Dw*mask==np.max(Dw*mask))
    print(ind_for_disp)
    return (pulsation_to_wavelength(np.max(Dw*mask)), n1_given[ind_for_disp],n2_given[ind_for_disp])

if __name__ == "__main__":
    target_central_wl = 600e-9
    target_puls = wavelength_to_pulsation(target_central_wl)
    epsilon = 10e-9