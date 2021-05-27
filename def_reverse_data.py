import math
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile
# First add all data to a list and then use the reversed 

def reverse(N_start, data_file):
    #REVERSE DATA
    #define counters
    N_total = N_start
    P_data = []

    for pressure in data_file:
        N_total = N_total + 1
        P_data.append(pressure)
    
    #define reversed data
    P_reversed = []

    for P in reversed(P_data) :
        P_reversed.append(P)

    return P_reversed, N_total