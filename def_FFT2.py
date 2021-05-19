import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile
from scipy import signal
import itertools

def ampli_phase_FFT(data_file, max_freqs, time):

    #define fft data
    fft_data = np.fft.fft(data_file)
    freq = np.fft.fftfreq(np.size(time), d = 0.001)

    #define lists for maxima intesities
    lst_magnitude = []
    lst_phase = []
    
    for f in max_freqs :

        index, = np.where(np.isclose(freq, f, atol = 0.001))

        magnitude = np.abs(fft_data[index[0]])
        phase = np.angle(fft_data[index[0]])
          
        lst_magnitude.append(magnitude)
        lst_phase.append(phase)

    print(f"\nFrequency: {max_freqs}\nMagnitude: {lst_magnitude} \nPhase: {lst_phase}\n ")

    return fft_data, freq, lst_magnitude, lst_phase