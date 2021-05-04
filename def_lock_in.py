import math
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile
import def_sine_wave
import def_cosine_wave
import def_butter_lowpass

def filter_freq(signal, frequency, max_freq, min_freq, alpha):

    #define amount of operations
    operations = len(signal)

    #define reference sine and cosine
    ref_sine = def_sine_wave.sine_wave(frequency, operations)
    ref_cosine = def_cosine_wave.cosine_wave(frequency, operations)

    #Multiply reference wave by signal and apply lowpass filter
    sin_mixed = np.multiply(signal, ref_sine)
    sin_mixed_lp = def_butter_lowpass.butter_lowpass_filter(sin_mixed, max_freq*0.8, alpha)

    cos_mixed = np.multiply(signal, ref_cosine)
    cos_mixed_lp = def_butter_lowpass.butter_lowpass_filter(cos_mixed, max_freq*0.8, alpha)

    #Define mean of sine and cosine wave
    mean_s = np.mean(sin_mixed_lp)
    mean_c = np.mean(cos_mixed_lp)

    #Define the magnitude and phase of the signal wave
    R = np.sqrt(np.add(np.square(cos_mixed_lp), np.square(sin_mixed_lp)))
    θ = np.arctan2(cos_mixed_lp, sin_mixed_lp)

    #calculate mean of magnitude and phase
    R_mean = np.mean(R)
    θ_mean = np.mean(θ)
    
    return R_mean, θ_mean  
    

        
    