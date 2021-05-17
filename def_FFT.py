import numpy as np
from scipy import signal
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
    max_freq_all = []
    intensity = maxima_freqs
    

    for x in itertools.repeat(0, (amount_freq + 2)) :
        del x

        max_index_int = np.argmax(intensity)
        max_intensity = intensity[max_index_int]

        intensity = np.where(intensity == max_intensity, 0, intensity)
        freq_index = freq_test_retrieval[max_index_int]
        max_freq_all.append(freq_index)

        if freq_index < (def_min_freq - 0.025):
            continue
        else :
            if len (max_freq) >= amount_freq:
                continue
            else: 
                max_freq.append(freq_index) 
        
    #sort frequencies
    max_freq.sort()
    
    return freq_test_retrieval, maxima_freqs, max_freq, f , Pxx_den




