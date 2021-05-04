import math
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile
from scipy import signal
import def_period

def complex_mod (θ, R_freq_P, R_freq_L,pip_rad, ):
    #INPUT:
    #   -phase difference θ_p-θ_l
    #   -oscillation amplitude ()
    #   -pipette radius
    #   -sample local radius
    beta1 = 2.0142
    beta3 = 2.1187
    
    
    
    #calculate the storage and loss modulus based on the extended
    #Zhou-Plaza method
    Estor = 3*pip_rad/beta1*R_freq_P/(1-(R/R_c)^beta3)/R_freq_L*cos(phi)
    Eloss = 3*pip_rad/beta1*R_freq_P/(1-(R/R_c)^beta3)/R_freq_L*sin(phi)
    tand = Eloss/Estor
    
    return Estor, Eloss, tand