import math
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile

def start_osci (N_wait, mean, standard_deviation , data_file, Fs):
    
    #count number of samples 
    N_samples_start = 0
    N_data_loop = N_wait + int(0.25*Fs)
    N_start = N_wait + int(0.25*Fs)

    #boolean to start and end new file
    booleanstart = False

    #define range of standard deviation
    mean_min = mean - (standard_deviation)
    mean_max = mean + (standard_deviation)

    #data for determination oscillation
    data_oscillation = data_file[N_data_loop:]
    
    for pressure in data_oscillation:

        #count how long signal is bigger that standard deviation. When the datapoints are bigger than waiting point + standard deviation for 0.25 seconds, N_start will be defined.
        # Again the boolean prevents the N_start to be overwritten and the mean_min en mean_max only change when datapoint is inbetween interval. 
        if booleanstart == False and (pressure > mean_max or pressure < mean_min) :
            N_samples_start = N_samples_start + 1

            #define start when signal is bigger for 0,025 seconds + boolean 
            if N_samples_start == int(0.25 * Fs):
                N_start = N_data_loop - int(0.25* Fs)
                booleanstart = True        
        else :
            N_samples_start = 0

            mean_min = pressure - 0.1*standard_deviation
            mean_max = pressure + 0.1*standard_deviation

        N_data_loop = N_data_loop + 1

    return (N_start)