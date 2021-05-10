import math
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile
import def_reverse_data

def stop_osci(N_wait_reversed, N_normal, standard_deviation, mean, data_file) :

    #define amount iterations
    N_samples_stop = 0
    N_now = N_wait_reversed 
    N_stop = N_wait_reversed

    #boolean to start and end new file
    booleanstop = False
    
    #define range of standard deviation
    mean_min = mean - (standard_deviation)
    mean_max = mean + (standard_deviation)
    
    #data for determination oscillation
    data_oscillation = data_file[N_normal:]

    length = len(data_file)
    N_end = length - N_normal
    stop_data = data_file[:N_end]

    #reversed data
    P_reversed, N_total = def_reverse_data.reverse(N_normal, stop_data)
    
    for pressure in P_reversed:

        #count how long signal is bigger that standard deviation
        if booleanstop == False and (pressure > mean_max or pressure < mean_min) :
            N_samples_stop = N_samples_stop + 1

            #define start when signal is bigger for 0,025 seconds + boolean 
            if N_samples_stop == 150:
                N_stop = N_now + 150
                booleanstop = True        
        else :
            N_samples_stop= 0
            mean_min = pressure - 0.1*standard_deviation
            mean_max = pressure + 0.1*standard_deviation

        N_now = N_now - 1
        
    return (N_stop)