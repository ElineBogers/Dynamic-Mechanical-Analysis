import math
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile
import os
import def_period
import def_reference
import def_FFT
import def_sine_wave
import def_cosine_wave
import def_lock_in
import def_butter_lowpass

#use Bayesian Methods for Hackers plotstyle. 
plt.style.use(["science", "no-latex", "grid", "vibrant"])

#decide which plots to show. when true it will be showed.
perio = False
creep = False

#define frequencies
frequencies = [0.05, 0.1, 0.25, 0.5, 1, 2, 6]

#to determine how many frequencies have to be found from local maxima
amount_freqs = len(frequencies)

#define phase and amplitude
amplitude = [50]*amount_freqs
phase = [0]*amount_freqs

#variable for how many times the data should be repeated
xN = 1

#how many cycle for the lowest frequency
Nlist=[1]*24        #, 2, 3, 4, 5, 6, 7, 8

#loop over graph properties
colors = ["Red", "Green", "Orange", "Purple", "Yellow", "Blue", "Brown", "Pink"]
x = 0

#define directory of files
directory = "Data_Eline"

#list for quality measurements
Ns = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10 , 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
areasIN = []
areasOUT = []
tick_label = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10" , "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"]

#deltasense data
for filename in os.listdir(str(directory)) :
    N = Nlist[x]
    
    if filename.endswith (".tdms") :
        tdms_file = TdmsFile.read(filename)
        group = tdms_file['Demodulated data']
        p_channel = group['pressure']

        #determine sampling frequency Fs from tdms file
        tdf = tdms_file["Raw data"]
        tdc = tdf.channels()
        Fs = int(tdc[0].properties["Effective_F_Sps"])

        #filter frequency
        p_channel = def_butter_lowpass.butter_lowpass_filter(p_channel, 8, Fs)


        #define functions
        time, oscillation, exta_oscillation = def_period.define_period(p_channel, xN, p_channel, perio, creep, Fs)

        #define amount of operations
        operations = len(oscillation)

        #import reference oscillation data input
        time1, signal = def_reference.reference_data(frequencies, amplitude, phase, operations, Fs)

        #define frequency by maxima of FFT
        f, maxima_f, max_freqs, f, Pxx_den  = def_FFT.maxima(oscillation, amount_freqs, xN, N, Fs)
        max_frequency = max(max_freqs)
        min_frequency = min(max_freqs)
    

        #Filter out single frequencies
        order = 2
        R, θ = def_lock_in.filter_freq(oscillation, max_freqs, order, Fs)

        #print(f" \n Frequencies: {max_freqs} \n Maginitude: {R} \n Phase: {θ} \n")

        #plot frequencies found
        t_gen, P_gen = def_reference.reference_data(max_freqs, R, θ, operations, Fs)
        
        #Plot reference input, signal and reference output
        plt.figure(2)
        #plt.plot(t_gen, P_gen, label = "Reference Output", color = "seagreen")
        plt.plot(time, oscillation, label = "Signal")
        plt.plot(time1, signal, label = "Reference Input")
        plt.legend(loc="upper right")
        plt.xlabel("time [s]")
        plt.ylabel("pressure [Pa]")

        plt.show()

        #starting value
        area_sig = 0
        area_input = 0
        area_output = 0
        area_in_out = 0

        #area difference
        for i in range(0,operations):
            area_sig = area_sig + abs(0.001 * ((signal[i]) + (signal[i - 1])) / 2)
            area_input = area_input + abs(0.001 * (((signal[i]) - (oscillation[i])) + ((signal[i - 1]) - (oscillation [i-1]))) / 2)
            area_output = area_output + abs(0.001 * (((P_gen[i]) - (oscillation[i])) + ((P_gen[i - 1])-(oscillation [i-1]))) / 2)
            area_in_out = area_in_out + abs((0.001 * ((P_gen[i]) + (P_gen[i - 1])) / 2) - (0.001 * ((signal[i]) + (signal[i - 1])) / 2))

        print(f"Area signal: {area_sig}, ∆ Area input: {area_input}, ∆ Area output: {area_output}, ∆ Area reference: {area_in_out}")

        #relative difference
        rela_IN = area_input/area_sig
        rela_OUT = area_output/area_sig

        #plot local maxima of FFT
        #plt.figure(2)
        #plt.loglog(f, Pxx_den)
        #plt.scatter(f, maxima_f)
        #plt.xlabel("frequency [Hz]")
        #plt.ylabel("Intensity")
        
        #append list to see area difference
        areasIN.append(rela_IN)
        areasOUT.append(rela_OUT)

        x = x + 1
        continue
    
    else :
        continue

    


plt.figure(3)
plt.bar(Ns, areasIN, tick_label = tick_label, width = 0.5, color = ["firebrick", 'royalblue'] )
plt.xlabel("N")
plt.ylabel("∆ Area")
plt.title("INPUT reference")

plt.figure(4)
plt.bar(Ns, areasOUT, tick_label = tick_label, width = 0.5, color = ["firebrick", 'royalblue'] )
plt.xlabel("N")
plt.ylabel("∆ Area")
plt.title("OUTPUT reference")

#plt.show()



