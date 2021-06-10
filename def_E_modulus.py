import math
import numpy as np
import matplotlib.pyplot as plt

#function to calculate elastic modulus, knowing the phase difference and amplitude of the DMA.

def elas_modulus(pip_rad, sam_rad, R_freq_P, R_freq_L, phi):
    E_storage = []
    E_loss = []
    Tand = []
    
    for i in range(0, len(R_freq_P)):
        #sample local radius
        beta1 = 2.0142
        beta3 = 2.1187

        #calculate the storage and loss modulus based on the extended
        #Zhou-Plaza method
        #assume 0.5 Poisson Ratio
        E_S = (3*pip_rad/beta1*R_freq_P[i]/(1-(pip_rad/sam_rad)**beta3)/R_freq_L[i]*math.cos(phi[i]))/10**3
        E_L = (3*pip_rad/beta1*R_freq_P[i]/(1-(pip_rad/sam_rad)**beta3)/R_freq_L[i]*math.sin(phi[i]))/10**3

        T = E_L/E_S

        E_storage.append(E_S)
        E_loss.append(E_L)
        Tand.append(T)

    return E_storage, E_loss, Tand