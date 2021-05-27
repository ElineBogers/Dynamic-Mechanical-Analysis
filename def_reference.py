import math
import numpy as np
import matplotlib.pyplot as plt

#create reference oscillation
def reference_data(frequency, amplitude, phase, operations, Fs) :
    #   input
    N_freq = len(frequency)

    #  pressure sensor takes 1 kHz samples (1000 samples/s)
    timestep = 1/Fs

    #   initialize command list
    signal = np.zeros(operations)
    time_signal = np.zeros(operations)

    #   do the actual work
    for i in range(0,operations):
        for n in range(0, N_freq):
               signal[i] = signal[i]+ amplitude[n]*math.sin((2*math.pi*frequency[n]*i*timestep) - phase[n])
        #define a time array for visualization purposes
        if i!=0: time_signal[i] = time_signal[i-1]+timestep

    signal = signal 
    return time_signal, signal


#
    #   calculate how many operations are needed
    #operations = int(min_test_duration/timestep)

