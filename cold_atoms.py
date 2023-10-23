import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, e, epsilon_0, c, m_u, k, m_e, alpha
from Reflectivity_our import *


def wavelength_to_pulsation(wl):
    return (2*np.pi*c)/(wl)

def pulsation_to_wavelength(om):
    return (2*np.pi*c)/(om)

def Susceptibility(x,y):
    return 3*((x)/(4*np.pi**2))*((2*y)+1j)/(1+4*(y)**2)

def Refraction_index(x,y):
    return np.real(np.sqrt(1+Susceptibility(x,y)))

def Absorption(x,y):
    return np.imag(np.sqrt(1+Susceptibility(x,y)))

def Reflection_array(density,n_air,d,lambda_atome,N):
    return np.array([reflectivity_new(n_air, Refraction_index(density,i), n_air, n_air, d*1e9, d*1e9, lambda_atome, N) for i in y])


if __name__ == "__main__":
    # Optical lattice
    wllaser = 800*1e-9
    d1 = wllaser/4
    N=1000
    n_air = 1

    # Atom parameters
    lambdas = np.array([770,780,800,830])
    x = np.linspace(1e-3,10,200)
    y = np.linspace(-5,5,200)

    for lambda_atome in lambdas:
        omega_atome = wavelength_to_pulsation(lambda_atome)
        #
        #------ Reflectivity colormap plot ------
        #
        #X,Y = np.meshgrid(x,y)
        #reflect = np.array([Reflection_array(density,n_air,d1,lambda_atome,N) for density in x]).T
        #plt.pcolor(Y,X,reflect, vmin=0, vmax=1)
        #plt.xlabel(r'$\Delta/\Gamma$')
        #plt.ylabel(r'$\rho_{at}\lambda^{3}$')
        #plt.title(r'$\lambda = $'+f'{lambda_atome} nm')
        #cbar = plt.colorbar()
        #cbar.set_label('Reflectivity')
        #plt.savefig(f'figures/reflectivity_colormap_lambda={lambda_atome}.pdf')
        #plt.close()

        # Now we fix the density and plot the reflectivity as a function of the detuning
        xfixed = 5
        
        
        # Plots of Reflectivity and absorption
        #
        reflect = Reflection_array(xfixed,n_air,d1,lambda_atome,N)
        absorption = Absorption(xfixed,y)
        plt.figure()
        fig, ax1 = plt.subplots()
        color = 'tab:red'
        ax1.set_xlabel(r'$\Delta/\Gamma$')
        ax1.set_ylabel('Reflectivity', color=color)
        ax1.plot(y,reflect,color=color)
        ax1.tick_params(axis='y', labelcolor=color)
        
        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
        
        color = 'tab:blue'
        ax2.set_ylabel(r'Optical extinction index $n$', color=color)  # we already handled the x-label with ax1
        ax2.plot(y,absorption, color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        
        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.title(r'$\lambda = $'+f'{lambda_atome} nm, '+r'$\rho_{at}\lambda^{3} = $'+f'{xfixed}')
        plt.savefig(f'figures/reflectivity_and_absorption_vs_detuning_lambda={lambda_atome}.pdf',bbox_inches='tight')
        plt.close()

        # Plot for refraction and absorption
        plt.figure()
        fig, ax1 = plt.subplots()
        color = 'tab:red'
        absorption = Absorption(xfixed,y)
        refraction = Refraction_index(xfixed,y)
        ax1.set_xlabel(r'$\Delta/\Gamma$')
        ax1.set_ylabel(r'Optical refractive index $n$', color=color)
        ax1.plot(y,refraction,color=color)
        ax1.tick_params(axis='y', labelcolor=color)
        
        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
        
        color = 'tab:blue'
        ax2.set_ylabel(r'Optical extinction index $\kappa$', color=color)  # we already handled the x-label with ax1
        ax2.plot(y,absorption, color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        
        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.title(r'$\lambda = $'+f'{lambda_atome} nm, '+r'$\rho_{at}\lambda^{3} = $'+f'{xfixed}')
        plt.savefig(f'figures/refraction_and_absorption_vs_detuning_lambda={lambda_atome}.pdf',bbox_inches='tight')
        plt.close()

