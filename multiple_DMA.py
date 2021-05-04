import math
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile
import os
import def_period
import def_FFT
import def_lock_in
import def_reference
import def_butter_lowpass

#to determine how many frequencies have to be found from local maxima
amount_freqs = 5

#define phase and amplitude
amplitude = [50]*amount_freqs
phase = [0]*amount_freqs

#variable for how many times the data should be repeated
xN = 1

#how many cycle for the lowest frequency
N = 1

#deltasense data
for filename in os.listdir(str("Data_Eline")) :
    if filename.endswith (".tdms") :
        tdms_file = TdmsFile.read(filename)
        group = tdms_file['Demodulated data']
        p_channel = group["pressure"]
        l_channel = group["aspirated length"] 
        
        #filter data with lowpass
        p_channel = def_butter_lowpass.butter_lowpass_filter(p_channel, 10)
        l_channel = def_butter_lowpass.butter_lowpass_filter(l_channel, 10)

        fig, ax1 = plt.subplots()
        ax1.plot(p_channel, label = "Pressure", color = 'royalblue')
        ax1.set_xlabel("operations")
        ax1.set_ylabel("Pressure [Pa]")
        plt.legend(loc="upper left")

        ax2 = ax1.twinx()
        ax2.plot(l_channel, label = "Aspirated length", color = "firebrick")
        ax2.set_ylabel("Displacement [nm]")
        plt.legend(loc="upper right")

        #Define oscillation (filter out preload)
        time, pressure, aspirated_length = def_period.define_period(p_channel, xN, l_channel)           

        #define amount of operations
        operations = len(pressure)
        
        #Defne displayed frequencies
        f, maxima_f, max_freqs, f, Pxx_den = def_FFT.maxima(pressure, amount_freqs, xN, N)
        max_frequency = max(max_freqs)
        min_frequency = min(max_freqs)

        #lists of magintude and phase
        R_P = []
        θ_P = []

        R_L = []
        θ_L = []

        θ_T = []

        #Filter out single frequencies
        for frequency in max_freqs:
            order = 2
            R_freq_P, θ_freq_P = def_lock_in.filter_freq(pressure, frequency, max_frequency, min_frequency, order)

            #append lists
            R_P.append(R_freq_P)
            θ_P.append(θ_freq_P)

            R_freq_L, θ_freq_L = def_lock_in.filter_freq(aspirated_length, frequency, max_frequency, min_frequency, order)

            #append lists
            R_L.append(R_freq_L)
            θ_L.append(θ_freq_L)

            θ = θ_freq_P - θ_freq_L

            θ_T.append(θ)

        
        

        print(f" \n Frequencies: {max_freqs} \n\n Maginitude Pressure: {R_P} \n Phase Pressure: {θ_P} \n\n Maginitude Aspirated Lenght: {R_L} \n Phase Aspirated Length: {θ_L} \n\n Phase difference: {θ_T}\n")

        #plot new amplitude and phase as a reference graph
        ref_P_time, ref_P = def_reference.reference_data(max_freqs, R_P, θ_P, operations)
        ref_L_time, ref_L = def_reference.reference_data(max_freqs, R_L, θ_L, operations)

        plt.show()

        plt.figure(3)
        plt.plot(time, pressure, label = "Pressure signal", color = 'royalblue')
        plt.plot(ref_P_time, ref_P, label = "Pressure reference", color = "firebrick")
        plt.xlabel("Time [s]")
        plt.ylabel("Pressure [Pa]")
        plt.legend(loc="upper left")

        plt.figure(4)
        plt.plot(time, aspirated_length, label = "Aspirated length signal", color = "royalblue")
        plt.plot(ref_L_time, ref_L, label = "Aspirated length reference", color = "firebrick")
        plt.xlabel("Time [s]")
        plt.ylabel("δOPL [nm]")
        plt.legend(loc="upper left")

        plt.show()

#
            #plt.figure(1)
            #plt.plot(R_freq, label = f"{frequency} Hz")
            #plt.legend(loc="upper left")
            #plt.title("Magnitude")
            
            #plt.figure(2)
            #plt.plot(θ_freq, label = f"{frequency} Hz")
            #plt.legend(loc="upper left")
            #plt.title("Phase") '''

        #plt.figure(1)
        #P_diff = np.diff(p_channel)
        #plt.plot(P_diff)
        #plt.show()


