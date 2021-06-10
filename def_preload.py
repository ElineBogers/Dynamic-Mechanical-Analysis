import math
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile

def define_wait(data_file, std_dev, Fs) :
    #define counters
    N_data = 5 * Fs
    N_equal = 0

    #define variables for pressure range
    P_prev1 = 0
    P_prev2 = 0

    #boolean for start wait
    booleanwait = False

    #define start wait
    N_wait = N_data

    #filter out preload (and time it takes to start run)
    start_data = data_file[N_data:]

    #loop over data to recognize wait. When pressure is still within a range of 7 times the standard deviation after 0.75 seconds, define start of wait.
    #Booleanwait is there to prevent to overwrite the time definition of the wait. Besides the boundaries of the interval are only changed when a datapoint is not within this range (P_prev1, P_prev2).
    for pressure in start_data :

        if P_prev1 < pressure < P_prev2 and booleanwait == False:
            #N_equal counts how many datapoints are within the interval consecutively.
            N_equal = N_equal + 1
            
            if N_equal == int(0.75*Fs):
                N_wait = N_data - int(0.75*Fs)
                booleanwait = True
        else:
            N_equal = 0

            P_prev1 = pressure - 8*std_dev
            P_prev2 = pressure + 8*std_dev

        N_data = N_data + 1

    #define end point for determination of standard deviation Interval of 0.5 seconds.
    N_wait_end = N_wait + int(0.5 *Fs)
    
    return (N_wait, N_wait_end)