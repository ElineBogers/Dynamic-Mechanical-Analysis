import math
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile
import os
from scipy import signal

from numpy.lib.twodim_base import triu_indices_from
import def_period
import def_FFT
import def_lock_in
import def_reference
import def_butter_lowpass
import def_E_modulus
import def_FFT2
from operator import truediv

#control which figure to show. When True it plots
perio = False
real = False
FFT = False
REF = False
Emod = False
Eav = True

#to determine how many frequencies have to be found from local maxima
amount_freqs = 7

#Sampling frquency
Fs = 4000

#define phase and amplitude
amplitude = [50]*amount_freqs
phase = [0]*amount_freqs

#variable for how many times the data should be repeated
xN = 1

#how many cycle for the lowest frequency
N = 2

#data foor generation fake data
frequencies = [5]
amplitude2 = [30]*amount_freqs
phase2 = [1/4*math.pi]*amount_freqs

#radius pipette and sample
rad_pip = 43e-6
rad_sam = 156e-6

#Counting amount of data files
amount = 0
E_storage_tot = [0] * amount_freqs
E_loss_tot = [0] * amount_freqs
E_storage_tot_FFT = [0] * amount_freqs
E_loss_tot_FFT = [0] * amount_freqs
lists_loss = {}
lists_storage = {}
lists_loss_FFT = {}
lists_storage_FFT = {}

for i in range (0,amount_freqs):
    lists_loss["f_loss"+str(i+1)] = []
    lists_storage["f_storage"+str(i+1)] = []
    lists_loss_FFT["f_loss_FFT"+str(i+1)] = []
    lists_storage_FFT["f_storage_FFT"+str(i+1)] = []

#deltasense data
for filename in os.listdir(str("Data_Eline")) :
    if filename.endswith (".tdms") :
        tdms_file = TdmsFile.read(filename)
        group = tdms_file['Demodulated data']
        p_channel = group["pressure"]
        l_channel = group["aspirated length"] 
        
        #filter data with lowpass
        p_channel = def_butter_lowpass.butter_lowpass_filter(p_channel, 10, Fs)
        l_channel = def_butter_lowpass.butter_lowpass_filter(l_channel, 10, Fs)

        #Define oscillation (filter out preload)
        time, pressure, aspirated_length = def_period.define_period(p_channel, xN, l_channel, perio, Fs)           
        
        #define amount of operations
        operations = len(pressure)
        
        #generation of fake data to see how program works. 
        #time, pressure = def_reference.reference_data(frequencies, amplitude, phase, operations)
        #Stime, aspirated_length = def_reference.reference_data(frequencies, amplitude2, phase2, operations)

        #Defne displayed frequencies
        f, maxima_f, max_freqs, f, Pxx_den = def_FFT.maxima(pressure, amount_freqs, xN, N, Fs)
        max_frequency = max(max_freqs)
        min_frequency = min(max_freqs)
        
        f_L, maxima_f_L, max_freqs_L, f_L, Pxx_den_L = def_FFT.maxima(aspirated_length, amount_freqs, xN, N, Fs)
        max_frequency = max(max_freqs)
        min_frequency = min(max_freqs)

        #new fft
        fft_data_P, freq_P, magnitude2_P, phase2_P = def_FFT2.ampli_phase_FFT(pressure, max_freqs, time, Fs)
        fft_data_L, freq_L, magnitude2_L, phase2_L = def_FFT2.ampli_phase_FFT(aspirated_length, max_freqs, time, Fs)
        

        #lists of magintude and phase
        R_P = []
        θ_P = []

        R_L = []
        θ_L = []

        θ_T = []
        θ_T_FFT = []

        E_storage = []
        E_loss = []
        Tand = []

        E_storage_FFT = []
        E_loss_FFT = []

        E_ratio = []
        E_ratio_FFT = []

        x = 0

        #Filter out single frequencies
        for frequency in max_freqs:
            order = 2
            R_freq_P, θ_freq_P = def_lock_in.filter_freq(pressure, frequency, max_frequency, order, Fs)

            #append lists
            R_P.append(R_freq_P)
            θ_P.append(θ_freq_P)

            R_freq_L, θ_freq_L = def_lock_in.filter_freq(aspirated_length, frequency, max_frequency, order, Fs)
            R_L_freq = R_freq_L / (1* math.exp(9))

            #append lists
            R_L.append(R_freq_L)
            θ_L.append(θ_freq_L)

            θ = (θ_freq_P - θ_freq_L)
            θ_T.append(θ)

            E_S, E_L, T = def_E_modulus.elas_modulus(rad_pip, rad_sam, R_freq_P, R_L_freq, θ)

            E_storage.append(E_S)
            E_loss.append(E_L)
            Tand.append(T)   

            magni_FFT_P = magnitude2_P[x]
            magni_FFT_L = magnitude2_L[x] / (1* math.exp(9))
            phase_FFT_P = phase2_P[x]
            phase_FFT_L = phase2_L[x]

            θ_FFT = (phase_FFT_P - phase_FFT_L)
            θ_T_FFT.append(θ_FFT)
            
            E_S_FFT, E_L_FFT, T_FFT = def_E_modulus.elas_modulus(rad_pip, rad_sam, magni_FFT_P, magni_FFT_L, θ_FFT)  #θ_FFT

            E_storage_FFT.append(E_S_FFT)
            E_loss_FFT.append(E_L_FFT)

            E_div = E_L/E_S
            E_div_FFT = E_L_FFT/E_S_FFT

            E_ratio.append(E_div)
            E_ratio_FFT.append(E_div_FFT)
            x = x + 1

        #plot new amplitude and phase as a reference graph
        ref_P_time, ref_P = def_reference.reference_data(max_freqs, R_P, θ_P, operations, Fs)
        ref_L_time, ref_L = def_reference.reference_data(max_freqs, R_L, θ_L, operations, Fs)

        ref_FFT_time_P, ref_FFT_P = def_reference.reference_data(max_freqs, magnitude2_P, phase2_P, operations, Fs)
        ref_FFT_time_L, ref_FFT_L = def_reference.reference_data(max_freqs, magnitude2_L, phase2_L, operations, Fs)

        print(f"\nFrequencies pressure: {max_freqs}\n\nFrequencies aspirated length: {max_freqs_L}  \n\nMaginitude Pressure: {R_P} \nPhase Pressure: {θ_P} \n\nMaginitude Aspirated Lenght: {R_L} \nPhase Aspirated Length: {θ_L} \n\nPhase difference: {θ_T}\n \nE': {E_storage}\nE'': {E_loss} \nE''/E' = {E_ratio} \n")
        print(f"\nMagnitude FFT Pressure: {magnitude2_P} \nPhase FFT pressure: {phase2_P}\n\nMagnitude FFT Aspirated Length: {magnitude2_L} \nPhase FFT Aspirated Length: {phase2_L}\n\nPhase difference FFT: {θ_T_FFT} \n\nE' FFT: {E_storage_FFT} \nE'' FFT: {E_loss_FFT}  \nE''/E' FFT: {E_ratio_FFT} \n")

        if real == True:

            fig4, ax1 = plt.subplots()
            ax1.plot(p_channel, label = "Pressure", color = 'royalblue')
            ax1.set_xlabel("operations")
            ax1.set_ylabel("Pressure [Pa]")
            plt.legend(loc="upper left")
            ax2 = ax1.twinx()
            ax2.plot(l_channel, label = "Aspirated length", color = "firebrick")
            ax2.set_ylabel("Displacement [nm]")
            plt.legend(loc="upper right")

            fig5, ax1 = plt.subplots()
            ax1.plot(l_channel, -p_channel)
            ax1.set_xlabel("Displacement [nm]")
            ax1.set_ylabel("Pressure [Pa]")

        if FFT == True:

            plt.figure(6)
            plt.loglog(f, Pxx_den)
            plt.title("FFT pressure")
            plt.xlabel("Frequency [Hz]")
            plt.ylabel("Intensity")

            plt.figure(7)
            plt.loglog(f_L, Pxx_den_L)
            plt.title("FFT aspirated length")
            plt.xlabel("Frequency [Hz]")
            plt.ylabel("Intensity")

        if REF == True:
            plt.figure(8)
            plt.plot(time, pressure, label = "Pressure signal", color = 'royalblue')
            plt.plot(ref_P_time, ref_P, label = "Pressure reference", color = "firebrick")
            plt.xlabel("Time [s]")
            plt.ylabel("Pressure [Pa]")
            plt.legend(loc="upper left")

            plt.figure(9)
            plt.plot(time, aspirated_length, label = "Aspirated length signal", color = "royalblue")
            plt.plot(ref_L_time, ref_L, label = "Aspirated length reference", color = "firebrick")
            plt.xlabel("Time [s]")
            plt.ylabel("δOPL [nm]")
            plt.legend(loc="upper left")

            plt.figure(10)
            plt.plot(time, pressure, label = "Pressure signal FFT", color = 'royalblue')
            plt.plot(ref_FFT_time_P, ref_FFT_P, label = "Pressure reference FFT", color = "firebrick")
            plt.xlabel("Time [s]")
            plt.ylabel("Pressure [Pa]")
            plt.legend(loc="upper left")

            plt.figure(11)
            plt.plot(time, aspirated_length, label = "Aspirated length signal FFT", color = "royalblue")
            plt.plot(ref_FFT_time_L, ref_FFT_L, label = "Aspirated length reference FFT", color = "firebrick")
            plt.xlabel("Time [s]")
            plt.ylabel("δOPL [nm]")
            plt.legend(loc="upper left")

            fig12, ax1 = plt.subplots()
            ax1.plot(aspirated_length, pressure, label = "Signal", color = 'royalblue')
            ax1.set_xlabel("Displacement [nm]")
            ax1.set_ylabel("Pressure [Pa]")
            plt.legend(loc="upper left")
            ax2 = ax1.twinx()
            ax2.plot(ref_FFT_L, ref_FFT_P, label = "Reference FFT", color = "firebrick")
            ax2.plot(ref_L, ref_P, label = "Reference Lock-In", color = "darkseagreen")
            ax2.set_ylabel("Pressure [Pa]")
            plt.legend(loc="upper right")
            
            

        if Emod == True:
            fig13, ax1 = plt.subplots()
            ax1.loglog(max_freqs, E_storage, label = "E'", color = "firebrick")
            ax1.scatter(max_freqs, E_storage, color = "firebrick")
            ax1.loglog(max_freqs, E_loss, label = "E''", color = "slateblue")
            ax1.scatter(max_freqs, E_loss, color = "slateblue")
            ax1.set_xlabel("Frequency [Hz]")
            ax1.set_ylabel("E', E'' [kPa]")
            plt.legend(loc="upper left")
            
            ax2 = ax1.twinx()
            ax2.plot(max_freqs, E_ratio, label = "tan δ", color = "goldenrod")
            ax2.scatter(max_freqs, E_ratio, color = "goldenrod")
            ax2.set_ylabel("tan δ")
            plt.legend(loc="upper right")
            plt.title("Lock-In determination")
           
            fig14, ax1 = plt.subplots()
            ax1.loglog(max_freqs, E_storage_FFT, label = "E' ", color = "firebrick")
            ax1.scatter(max_freqs, E_storage_FFT, color = "firebrick")
            ax1.loglog(max_freqs, E_loss_FFT, label = "E''", color = "slateblue")
            ax1.scatter(max_freqs, E_loss_FFT, color = "slateblue")
            ax1.set_xlabel("Frequency [Hz]")
            ax1.set_ylabel("E', E'' [kPa]")
            plt.legend(loc="upper left")
            
            ax2 = ax1.twinx()
            ax2.plot(max_freqs, E_ratio_FFT, label = "tan δ", color = "goldenrod")
            ax2.scatter(max_freqs, E_ratio_FFT, color = "goldenrod")
            ax2.set_ylabel("tan δ")
            plt.legend(loc="upper right")
            plt.title("FFT determination")


        plt.show()

       
    for i in range(0, amount_freqs):
        lists_loss["f_loss"+str(i+1)].append(E_loss[i])
        lists_storage["f_storage"+str(i+1)].append(E_storage[i])

        lists_loss_FFT["f_loss_FFT"+str(i+1)].append(E_loss_FFT[i])
        lists_storage_FFT["f_storage_FFT"+str(i+1)].append(E_storage_FFT[i])
    
    E_storage_tot = np.add(E_storage_tot, E_storage)
    E_loss_tot = np.add(E_loss_tot, E_loss)
    E_storage_tot_FFT = np.add(E_storage_tot_FFT, E_storage_FFT)
    E_loss_tot_FFT = np.add(E_loss_tot_FFT, E_loss_FFT)

    amount = amount + 1
    
E_storage_mean = E_storage_tot / amount
E_loss_mean = E_loss_tot / amount

E_storage_mean_FFT = E_storage_tot_FFT / amount
E_loss_mean_FFT = E_loss_tot_FFT / amount

Eloss_std_dev = []
Estorage_std_dev = []
Eloss_std_dev_FFT = []
Estorage_std_dev_FFT = []

for freq in range(0, amount_freqs):
    std_dev_loss = lists_loss["f_loss"+str(freq + 1)]
    std_dev_storage = lists_storage["f_storage"+str(freq + 1)]
    std_dev_loss_FFT = lists_loss_FFT["f_loss_FFT"+str(freq + 1)]
    std_dev_storage_FFT = lists_storage_FFT["f_storage_FFT"+str(freq + 1)]

    lst_std_dev_loss = []
    lst_std_dev_storage = []
    lst_std_dev_loss_FFT = []
    lst_std_dev_storage_FFT = []

    for x in range(0, amount):
        std_dev_L = (std_dev_loss[x] - E_loss_mean[freq]) ** 2
        std_dev_S = (std_dev_storage[x] - E_storage_mean[freq]) ** 2
        std_dev_L_FFT = (std_dev_loss_FFT[x] - E_loss_mean_FFT[freq]) ** 2
        std_dev_S_FFT = (std_dev_storage_FFT[x] - E_storage_mean_FFT[freq]) ** 2

        lst_std_dev_loss.append(std_dev_L)
        lst_std_dev_storage.append(std_dev_S)
        lst_std_dev_loss_FFT.append(std_dev_L_FFT)
        lst_std_dev_storage_FFT.append(std_dev_S_FFT)

    std_dev_L_fin = math.sqrt(sum(lst_std_dev_loss) / amount)
    std_dev_S_fin = math.sqrt(sum(lst_std_dev_storage) / amount)
    std_dev_L_fin_FFT = math.sqrt(sum(lst_std_dev_loss_FFT) / amount)
    std_dev_S_fin_FFT = math.sqrt(sum(lst_std_dev_storage_FFT) / amount)

    Eloss_std_dev.append(std_dev_L_fin)
    Estorage_std_dev.append(std_dev_S_fin)
    Eloss_std_dev_FFT.append(std_dev_L_fin_FFT)
    Estorage_std_dev_FFT.append(std_dev_S_fin_FFT)

print(f"\nEloss variance: {Eloss_std_dev}  \nEstorage variance: {Estorage_std_dev} \n\nEloss standard deviation FFT: {Eloss_std_dev_FFT}  \nEstorage standard deviation FFT: {Estorage_std_dev_FFT}")

if Eav == True:

    plt.figure(15)
    plt.loglog(max_freqs, E_storage_tot, label = "E'", color = "firebrick")
    plt.scatter(max_freqs, E_storage_tot, color = "firebrick")
    plt.errorbar(max_freqs, E_storage_tot, yerr = Estorage_std_dev, color = "firebrick")
    plt.loglog(max_freqs, E_loss_tot, label = "E''", color = "royalblue")
    plt.scatter(max_freqs, E_loss_tot, color = "royalblue")
    plt.errorbar(max_freqs, E_loss_tot, yerr = Eloss_std_dev, color = "royalblue")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("E', E'' [kPa]")
    plt.legend(loc="upper left")
    plt.title(f"Elastic modulus on N = {amount} determined with Lock-In")

    plt.figure(16)
    plt.loglog(max_freqs, E_storage_tot_FFT, label = "E' FFT", color = "firebrick")
    plt.scatter(max_freqs, E_storage_tot_FFT, color = "firebrick")
    plt.errorbar(max_freqs, E_storage_tot_FFT, yerr = Estorage_std_dev_FFT, color = "firebrick")
    plt.loglog(max_freqs, E_loss_tot_FFT, label = "E'' FFT", color = "royalblue")
    plt.scatter(max_freqs, E_loss_tot_FFT, color = "royalblue")
    plt.errorbar(max_freqs, E_loss_tot_FFT, yerr = Eloss_std_dev_FFT, color = "royalblue")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("E', E'' [kPa]")
    plt.legend(loc="upper left")
    plt.title(f"Elastic modulus on N = {amount} determined with FFT")

    corr = signal.correlate(ref_L,aspirated_length, mode = "full")
    lags = signal.correlation_lags(len(ref_L), len(aspirated_length), mode = "full")
    corr /= np.max(corr)
    plt.figure(20)
    plt.plot(lags,corr)

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


