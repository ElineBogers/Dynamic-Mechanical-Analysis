import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile
from scipy import signal
import x_def_bandpass
import itertools
import math

def ampli_phase_FFT(data_file, max_freqs, time, Fs):

    #define fft data
    fft_data = np.fft.fft(data_file)
    freq = np.fft.fftfreq(np.size(time), d = 1/Fs)

    #define lists for maxima intesities
    lst_magnitude = []
    lst_phase = []

    N = len(data_file)
    
    #define magnitude and phase for each frequency.
    for f in max_freqs :
        
        fft_data_window = np.fft.fft(data_file*signal.windows.gaussian(len(data_file), std = 10, sym = True))

        index, = np.where(np.isclose(freq, f, atol = 0.001))
        magnitude = np.abs(fft_data[index[0]]) / N * 2

        phase = np.angle(fft_data_window[index[0]])
          
        lst_magnitude.append(magnitude)
        lst_phase.append(phase)         #+ 1/2*math.pi

    return fft_data, freq, lst_magnitude, lst_phase


    #N = 1
        #periods_freq = int(1/f * N *1000)
        #operations_tot = len(data_file)

        #while periods_freq > operations_tot:
         #   N = N - 1
          #  periods_freq = int(1/f * N *1000)

        #else:
           # determine_data = data_file[0:periods_freq]
           # time = time[0:periods_freq]
       # print(periods_freq)