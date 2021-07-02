import math
import numpy as np
import matplotlib.pyplot as plt
import def_butter_bandpass
import def_butter_lowpass
import sympy

def lock_in_custom(data, frequencies, order, Fs):

    #list for amplitude and phase of iets frequency
    list_R = []
    list_θ = []

    for frequency in frequencies:

        #create reference wave by filtering out single frequency
        data_single_sin = def_butter_bandpass.butter_bandpass_filter(data, (frequency-0.01), (frequency+0.01), Fs, 1)
        
        #create cosine of reference wave by calculating the derivative.
        data_single_cos = np.diff(data_single_sin) / (1/Fs)
        data_single_cos = np.append(data_single_cos, 0)

        #Multiply reference waves by signal and apply lowpass filter
        sin_mixed = np.multiply(data, data_single_sin)
        sin_mixed_lp = def_butter_lowpass.butter_lowpass_filter(sin_mixed, frequency*1, Fs, order)

        cos_mixed = np.multiply(data, data_single_cos)
        cos_mixed_lp = def_butter_lowpass.butter_lowpass_filter(cos_mixed, frequency*1, Fs, order)

        #Define the magnitude and phase of the signal wave
        R = np.sqrt(np.add(np.square(cos_mixed_lp), np.square(sin_mixed_lp)))
        θ = np.arctan2(cos_mixed_lp, sin_mixed_lp)

        #calculate mean of magnitude and phase
        R_mean = np.mean(R)
        θ_mean = np.mean(θ)

        list_R.append(R_mean)
        list_θ.append(θ_mean)

    
    return list_R, list_θ 
    