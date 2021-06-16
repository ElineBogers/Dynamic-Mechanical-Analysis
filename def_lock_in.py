import math
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile
import def_sine_wave
import def_cosine_wave
import def_butter_lowpass

def filter_freq(signal, frequencies, order, Fs):

    list_R = []
    list_θ = []
    min_freq = min(frequencies)

    for frequency in frequencies:

        #define period where to extract the information from. N gives the amount of periods that are included per frequency.
        N = 100
        periods_freq = int(1/frequency * N *Fs)
        operations_tot = len(signal)

        #keep on removing periods untill the length of the amount of periods is equal or smaller than 
        while periods_freq > operations_tot:
            N = N - 1
            periods_freq = int(1/frequency * N *Fs)

        else:
            determine_data = signal[0:periods_freq]

        #define amount of operations
        operations = len(determine_data)

        #define reference sine and cosine
        ref_sine = def_sine_wave.sine_wave(frequency, operations)
        ref_cosine = def_cosine_wave.cosine_wave(frequency, operations)

        #Multiply reference wave by signal and apply lowpass filter
        sin_mixed = np.multiply(determine_data, ref_sine)
        sin_mixed_lp = def_butter_lowpass.butter_lowpass_filter(sin_mixed, min_freq, Fs, order)

        cos_mixed = np.multiply(determine_data, ref_cosine)
        cos_mixed_lp = def_butter_lowpass.butter_lowpass_filter(cos_mixed, min_freq, Fs, order)

        #Define mean of sine and cosine wave
        mean_s = np.mean(sin_mixed_lp)
        mean_c = np.mean(cos_mixed_lp)

        #Define the magnitude and phase of the signal wave
        R = np.sqrt(np.add(np.square(cos_mixed_lp), np.square(sin_mixed_lp)))
        θ = np.arctan2(cos_mixed_lp, sin_mixed_lp)

        #calculate mean of magnitude and phase
        R_mean = np.mean(R)
        θ_mean = np.mean(θ)

        list_R.append(R_mean)
        list_θ.append(θ_mean)
    
    return list_R, list_θ 
    

        
    