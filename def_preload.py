import math
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile

def define_wait(data_file, std_dev) :
    #define counters
    N_data = 5000
    N_equal = 0

    #define variables for pressure range
    P_prev1 = 0
    P_prev2 = 0

    #boolean for start wait
    booleanwait = False

    #define start wait
    N_wait = 0

    #filter out preload (and time it takes to start run)
    start_data = data_file[5000:]

    #loop over data to recognize wait
    for pressure in start_data :

        if P_prev1 < pressure < P_prev2 and booleanwait == False:
            N_equal = N_equal + 1
            
            if N_equal == 750:
                N_wait = N_data - 750
                booleanwait = True
        else:
            N_equal = 0

            P_prev1 = pressure - 7*std_dev
            P_prev2 = pressure + 7*std_dev

        N_data = N_data + 1

    #define end point for determination of standard deviation
    N_wait_end = N_wait + 500
    #print (N_wait)
    return (N_wait, N_wait_end)