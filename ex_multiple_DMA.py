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
import def_E_modulus
import def_FFT2

#to determine how many frequencies have to be found from local maxima
amount_freqs = 6

#define phase and amplitude
amplitude = [50]*amount_freqs
phase = [0]*amount_freqs

#variable for how many times the data should be repeated
xN = 1

#how many cycle for the lowest frequency
N = 2

#radius pipette and sample
rad_pip = 86e-6
rad_sam = 350e-6

#Counting amount of data files
amount = 0
E_storage_tot = [0] * amount_freqs
E_loss_tot = [0] * amount_freqs
lists_loss = {}
lists_storage = {}

for i in range (0,amount_freqs):
    lists_loss["f_loss"+str(i+1)] = []
    lists_storage["f_storage"+str(i+1)] = []

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

        #Define oscillation (filter out preload)
        time, pressure, aspirated_length = def_period.define_period(p_channel, xN, l_channel)           

        #define amount of operations
        operations = len(pressure)
        
        #Defne displayed frequencies
        f, maxima_f, max_freqs, f, Pxx_den = def_FFT.maxima(pressure, amount_freqs, xN, N)
        max_frequency = max(max_freqs)
        min_frequency = min(max_freqs)

        f_L, maxima_f_L, max_freqs_L, f_L, Pxx_den_L = def_FFT.maxima(aspirated_length, amount_freqs, xN, N)
        max_frequency = max(max_freqs)
        min_frequency = min(max_freqs)

        #new fft
        fft_data_P, freq_P, magnitude2_P, phase2_P = def_FFT2.ampli_phase_FFT(pressure, max_freqs, time)
        fft_data_L, freq_L, magnitude2_L, phase2_L = def_FFT2.ampli_phase_FFT(aspirated_length, max_freqs_L, time)


        #lists of magintude and phase
        R_P = []
        θ_P = []

        R_L = []
        θ_L = []

        θ_T = []

        E_storage = []
        E_loss = []
        Tand = []

        E_storage_FFT = []
        E_loss_FFT = []
        x = 0

        #Filter out single frequencies
        for frequency in max_freqs:
            order = 2
            R_freq_P, θ_freq_P = def_lock_in.filter_freq(pressure, frequency, max_frequency, min_frequency, order)

            #append lists
            R_P.append(R_freq_P)
            θ_P.append(θ_freq_P)

            R_freq_L, θ_freq_L = def_lock_in.filter_freq(aspirated_length, frequency, max_frequency, min_frequency, order)
            R_L_freq = R_freq_L / (1* math.exp(9))

            #append lists
            R_L.append(R_freq_L)
            θ_L.append(θ_freq_L)

            θ = θ_freq_P - θ_freq_L

            θ_T.append(θ)

            E_S, E_L, T = def_E_modulus.elas_modulus(rad_pip, rad_sam, R_freq_P, R_L_freq, θ)

            E_storage.append(E_S)
            E_loss.append(E_L)
            Tand.append(T)   

            magni_FFT_P = magnitude2_P[x]
            magni_FFT_L = magnitude2_L[x]
            phase_FFT_P = phase2_P[x]
            phase_FFT_L = phase2_L[x]

            θ_FFT = phase_FFT_P - phase_FFT_L
            
            E_S_FFT, E_L_FFT, T_FFT = def_E_modulus.elas_modulus(rad_pip, rad_sam, magni_FFT_P, magni_FFT_L, θ_FFT)

            E_storage_FFT.append(E_S_FFT)
            E_loss_FFT.append(E_L_FFT)

            x = x + 1
        
        print(f" \n Frequencies pressure: {max_freqs}\n\n Frequencies aspirated length: {max_freqs_L}  \n\n Maginitude Pressure: {R_P} \n Phase Pressure: {θ_P} \n\n Maginitude Aspirated Lenght: {R_L} \n Phase Aspirated Length: {θ_L} \n\n Phase difference: {θ_T}\n \n E': {E_storage}\n E'': {E_loss} \n \n E' FFT: {E_storage_FFT} \n E'' FFT: {E_loss_FFT}")

        #plot new amplitude and phase as a reference graph
        ref_P_time, ref_P = def_reference.reference_data(max_freqs, R_P, θ_P, operations)
        ref_L_time, ref_L = def_reference.reference_data(max_freqs, R_L, θ_L, operations)

        ref_FFT_time_P, ref_FFT_P = def_reference.reference_data(max_freqs, magnitude2_P, phase2_P, operations)
        ref_FFT_time_L, ref_FFT_L = def_reference.reference_data(max_freqs, magnitude2_L, phase2_L, operations)

        fig, ax1 = plt.subplots()
        ax1.plot(p_channel, label = "Pressure", color = 'royalblue')
        ax1.set_xlabel("operations")
        ax1.set_ylabel("Pressure [Pa]")
        plt.legend(loc="upper left")

        ax2 = ax1.twinx()
        ax2.plot(l_channel, label = "Aspirated length", color = "firebrick")
        ax2.set_ylabel("Displacement [nm]")
        plt.legend(loc="upper right")

        plt.figure(2)
        plt.loglog(f, Pxx_den)
        plt.title("FFT pressure")
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Intensity")

        plt.figure(3)
        plt.loglog(f_L, Pxx_den_L)
        plt.title("FFT aspirated length")
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Intensity")

        plt.figure(4)
        plt.plot(time, pressure, label = "Pressure signal", color = 'royalblue')
        plt.plot(ref_P_time, ref_P, label = "Pressure reference", color = "firebrick")
        plt.xlabel("Time [s]")
        plt.ylabel("Pressure [Pa]")
        plt.legend(loc="upper left")

        plt.figure(5)
        plt.plot(time, aspirated_length, label = "Aspirated length signal", color = "royalblue")
        plt.plot(ref_L_time, ref_L, label = "Aspirated length reference", color = "firebrick")
        plt.xlabel("Time [s]")
        plt.ylabel("δOPL [nm]")
        plt.legend(loc="upper left")

        plt.figure(6)
        plt.loglog(max_freqs, E_storage, label = "E'", color = "orange")
        plt.scatter(max_freqs, E_storage)
        plt.loglog(max_freqs, E_loss, label = "E''", color = "slateblue")
        plt.scatter(max_freqs, E_loss)
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("E', E'' [kPa]")
        plt.legend(loc="upper left")


        plt.figure(7)
        plt.plot(time, pressure, label = "Pressure signal", color = 'royalblue')
        plt.plot(ref_FFT_time_P, ref_FFT_P, label = "Pressure reference", color = "firebrick")
        plt.xlabel("Time [s]")
        plt.ylabel("Pressure [Pa]")
        plt.legend(loc="upper left")

        plt.figure(8)
        plt.plot(time, aspirated_length, label = "Aspirated length signal FFT", color = "royalblue")
        plt.plot(ref_FFT_time_L, ref_FFT_L, label = "Aspirated length reference FFT", color = "firebrick")
        plt.xlabel("Time [s]")
        plt.ylabel("δOPL [nm]")
        plt.legend(loc="upper left")

        plt.figure(9)
        plt.loglog(max_freqs, E_storage_FFT, label = "E' FFT", color = "orange")
        plt.scatter(max_freqs, E_storage_FFT)
        plt.loglog(max_freqs, E_loss_FFT, label = "E'' FFT", color = "slateblue")
        plt.scatter(max_freqs, E_loss_FFT)
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("E', E'' [kPa]")
        plt.legend(loc="upper left")

        plt.show()

       
    for i in range(0, amount_freqs):
        lists_loss["f_loss"+str(i+1)].append(E_loss[i])
        lists_storage["f_storage"+str(i+1)].append(E_storage[i])
    
    E_storage_tot = np.add(E_storage_tot, E_storage)
    E_loss_tot = np.add(E_loss_tot, E_loss)
    amount = amount + 1
    
E_storage_mean = E_storage_tot / amount
E_loss_mean = E_loss_tot / amount

Eloss_variance = []
Estorage_variance = []

for freq in range(0, amount_freqs):
    std_dev_loss = lists_loss["f_loss"+str(freq + 1)]
    std_dev_storage = lists_storage["f_storage"+str(freq + 1)]

    variances_loss = []
    variances_storage = []

    for x in range(0, amount):
        std_dev_L = (std_dev_loss[x] - E_loss_mean[freq]) ** 2
        std_dev_S = (std_dev_storage[x] - E_storage_mean[freq]) ** 2

        variances_loss.append(std_dev_L)
        variances_storage.append(std_dev_S)

    variance_L = math.sqrt(sum(variances_loss) / amount)
    variance_S = math.sqrt(sum(variances_loss) / amount)

    Eloss_variance.append(variance_L)
    Estorage_variance.append(variance_S)

plt.loglog(max_freqs, E_storage_tot, label = "E' FFT", color = "orange")
plt.scatter(max_freqs, E_storage_tot)
plt.errorbar(max_freqs, E_storage_tot, yerr = Estorage_variance)
plt.loglog(max_freqs, E_loss_tot, label = "E'' FFT", color = "slateblue")
plt.scatter(max_freqs, E_loss_tot)
plt.errorbar(max_freqs, E_loss_tot, yerr = Eloss_variance)
plt.xlabel("Frequency [Hz]")
plt.ylabel("E', E'' [kPa]")
plt.legend(loc="upper left")
plt.title(f"E modulus on N = {amount}")
plt.show()

print(f"Eloss variance: {Eloss_variance}  \nEstorage variance: {Estorage_variance}")

 

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


