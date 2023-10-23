# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 16:00:55 2023

@author: jbaza
"""
import numpy as np
import matplotlib.pyplot as plt 

def matrix (n, d, lm):
    
   
    M = np.array([[np.cos((2*np.pi*n)*d/lm), -(1j/n)*np.sin((2*np.pi*n)*d/lm)],[(-1j*n)*np.sin((2*np.pi*n)*d/lm), np.cos((2*np.pi*n)*d/lm)]])
    return M

def reflectivity(n_0, n_1, n_2, n_s, d_1, d_2, lm, N):
    
    M = np.dot(matrix(n_1, d_1, lm), matrix(n_2, d_2, lm))
    result = np.linalg.matrix_power(M,N)
    reflect = pow(np.abs((result[0,0]*n_0+result[0,1]*n_0*n_s - result[1,0] - result[1,1])/(result[0,0]*n_0+result[0,1]*n_0*n_s + result[1,0] + result[1,1])),2)
    return reflect

def transmittivity(n_0, n_1, n_2, n_s, d_1, d_2, lm, N):
    
    M = np.dot(matrix(n_1, d_1, lm), matrix(n_2, d_2, lm))
    result = np.linalg.matrix_power(M,N)
    transmit = pow(np.abs(2*n_0/(result[0,0]*n_0+result[0,1]*n_0*n_s + result[1,0] + result[1,1])),2)
    return transmit


def reflectivity_new(n_0, n_1, n_2, n_s, d_1, d_2, lm, N):
    
    return 1-np.array(transmittivity(n_0, n_1, n_2, n_s, d_1, d_2, lm, N))

def function():
    X = np.linspace(400,900, 400)  
    #Y = [reflectivity(4.0,4.0,3.32,3.32,37.0,44.0,i, 20 ) for i in X]
    Z = [transmittivity(4.0,4.0,3.32,3.32,37.0,44.0,i, 20) for i in X]
    R = 1 - np.array(Z)
    plt.figure()
    plt.plot(X, R)
    #plt.plot(X, Z)
    # plt.plot(X,S)
    plt.show()
    plt.close()
   
  
if __name__ == '__main__':
    function()

    
    
    
    
        
        