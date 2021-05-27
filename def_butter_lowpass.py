import numpy as np
from scipy.signal import butter, lfilter, freqz
import matplotlib.pyplot as plt

#Fiter that is cutting of frequencies higher than the cutoff frequency.

def butter_lowpass(cutoff, Fs, order=5):
    nyq = 0.5 * Fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, Fs, order=5):
    b, a = butter_lowpass(cutoff, Fs, order=order)
    y = lfilter(b, a, data)
    return y