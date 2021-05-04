import math
import numpy as np
import matplotlib.pyplot as plt


def elas_modulus(pip_rad, sam_rad, R_freq_P, R_freq_L, phi):

    #sample local radius
    beta1 = 2.0142
    beta3 = 2.1187

    #calculate the storage and loss modulus based on the extended
    #Zhou-Plaza method
    #assume 0.5 Poisson Ratio
    Estor = 3*pip_rad/beta1*R_freq_P/(1-(pip_rad/sam_rad)**beta3)/R_freq_L*math.cos(phi)
    Eloss = 3*pip_rad/beta1*R_freq_P/(1-(pip_rad/sam_rad)**beta3)/R_freq_L*math.sin(phi)

    tand = Eloss/Estor

    return Estor, Eloss, tand