import math
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile

def define_wait(data_file, std_dev) :
    #define counters
    N_data = 4000
    N_equal = 0

    #define variables for pressure range
    P_prev1 = 0
    P_prev2 = 0

    #boolean for start wait
    booleanwait = False

    #define start wait
    N_wait = 0

    #filter out preload (and time it takes to start run)
    start_data = data_file[4000:]

    #loop over data to recognize wait
    for pressure in start_data :

        if P_prev1 < pressure < P_prev2 and booleanwait == False:
            N_equal = N_equal + 1
            
            if N_equal == 150:
                N_wait = N_data - 150
                booleanwait = True
        else:
            N_equal = 0

            P_prev1 = pressure - std_dev
            P_prev2 = pressure + std_dev

        N_data = N_data + 1

    #define end point for determination of standard deviation
    N_wait_end = N_wait + 1000
    #print (N_wait)
    return (N_wait, N_wait_end)