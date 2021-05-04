import math
import numpy as np
import matplotlib.pyplot as plt

def sine_wave(freq, operations):

    #define list for signal
    signal = np.zeros(operations)

    timestep = 1/1000

    #loop over operations to define signal
    for i in range(0,operations):
        #define signal array
        signal[i] = signal[i] + math.sin(2*math.pi*freq*i*timestep)


    return signal   

#    
    #define how long test takes
    #"min_test_duration = 1/min_freq * N_min  "

    #pressure sensor takes 1 kHz samples (1000 samples/s)
    #"timestep = 1/1000 "

    #calculate how many operations are needed
    #"operations = int(min_test_duration/timestep) "  
