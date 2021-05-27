import numpy as np
from scipy import signal
import itertools
import matplotlib.pyplot as plt

def maxima (data_file, amount_freq, xN, N, Fs):
    len_data = len(data_file)

    #for x in range(0,len_data):
        #data_file.append(0)

    f, Pxx_den = signal.periodogram(data_file, Fs)

    #define amount operations and minimal frequency
    #operations = len(data_file)
    def_min_freq = 1/(len_data * (1/Fs) / (N*xN))

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

        if freq_index < (def_min_freq - 0.01):
            continue
        else :
            if len (max_freq) >= amount_freq:
                continue
            else: 
                max_freq.append(freq_index) 
        
    #sort frequencies
    max_freq.sort()
    
    return freq_test_retrieval, maxima_freqs, max_freq, f , Pxx_den