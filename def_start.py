import math
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile

def start_osci (N_wait, mean, standard_deviation , data_file):
    
    #count number of samples 
    N_samples_start = 0
    N_data_loop = N_wait
    N_start = int()

    #boolean to start and end new file
    booleanstart = False

    #define range of standard deviation
    mean_min = mean - (standard_deviation)
    mean_max = mean + (standard_deviation)

    #data for determination oscillation
    data_oscillation = data_file[N_wait:]
    
    for pressure in data_oscillation:

        #count how long signal is bigger that standard deviation
        if booleanstart == False and (pressure > mean_max or pressure < mean_min) :
            N_samples_start = N_samples_start + 1

            #define start when signal is bigger for 0,025 seconds + boolean 
            if N_samples_start == 150:
                N_start = N_data_loop - 150
                booleanstart = True        
        else :
            N_samples_start = 0

            mean_min = pressure - 0.1*standard_deviation
            mean_max = pressure + 0.1*standard_deviation

        N_data_loop = N_data_loop + 1
    return (N_start)