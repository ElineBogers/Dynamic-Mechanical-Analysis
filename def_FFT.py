import math
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile
from scipy import signal
import def_period
from scipy.signal import find_peaks
import itertools

def maxima (data_file, amount_freq, xN, N):
    #define timestep (1kHz sampling freq)
    fs = 1000 

    f, Pxx_den = signal.periodogram(data_file, fs)

    #define amount operations and minimal frequency
    operations = len(data_file)
    def_min_freq = 1/(operations * (1/1000) / (N*xN))

    #find local maxima
    freq_test_retrieval = f[signal.argrelextrema(Pxx_den,np.greater)]
    maxima_freqs = Pxx_den[signal.argrelextrema(Pxx_den,np.greater)]

    #define lists for maxima intesities
    max_freq = [] * amount_freq
    intensity = maxima_freqs

    for x in itertools.repeat(0, (amount_freq + 2)) :
        del x

        max_index_int = np.argmax(intensity)
        max_intensity = intensity[max_index_int]

        intensity = np.where(intensity == max_intensity, 0, intensity)
        freq_index = freq_test_retrieval[max_index_int]

        if freq_index < (def_min_freq - 0.001):
            continue
        else :
            if len (max_freq) >= amount_freq:
                continue
            else: 
                max_freq.append(freq_index) 

    #define minimum frequency by time of oscillation
    min_freq = min(max_freq)
    max_freq = [def_min_freq if i == min_freq else i for i in max_freq]
    
    plt.figure(2)
    plt.loglog(f, Pxx_den)
    plt.title("FFT")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Intensity")

    
    return freq_test_retrieval, maxima_freqs, max_freq, f , Pxx_den




