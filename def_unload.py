import math
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile
import def_reverse_data

def define_unload(data_file, std_dev):

    #define variables for pressure range
    P_prev1 = 0
    P_prev2 = 0

    #boolean for start wait
    booleanwait = False

    #define start wait
    N_wait_reversed = 0

    #filter out preload (and time it takes to start run)
    N_start = 4000
    length = len(data_file)
    N_end = length - N_start
    start_data = data_file[:N_end]

    #reversed data
    P_reversed, N_total = def_reverse_data.reverse(N_start, start_data)
    
    #define counters
    N_equal = 0
    N_now = N_total - N_start

    #loop over data to recognize wait
    for pressure in P_reversed :

        if P_prev1 < pressure < P_prev2 and booleanwait == False:
            N_equal = N_equal + 1
            
            if N_equal == 500:
                N_wait_reversed = N_now + 300
                booleanwait = True
        else:
            N_equal = 0

            P_prev1 = pressure - 5*std_dev
            P_prev2 = pressure + 5*std_dev
        N_now = N_now - 1

    N_normal = N_total - N_wait_reversed
    N_wait_reversed_start = N_wait_reversed - 500

    return N_wait_reversed, N_normal, N_wait_reversed_start
    

    
