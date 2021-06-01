import math
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile
import os
from numpy.lib.function_base import hamming
from scipy import signal

from numpy.lib.twodim_base import triu_indices_from
from scipy.signal.windows.windows import hann
import def_period
import def_FFT
import def_lock_in
import def_reference
import def_butter_lowpass
import def_E_modulus
import def_FFT2
import x_def_bandpass
from operator import truediv

#use Bayesian Methods for Hackers plotstyle. 
plt.style.use('bmh')

#control which figure to show. When True plot is shown
perio = True           #shows whole measurement with segmentation of different stages
real = True            #shows both aspirated length and pressure sample results
FFT = True             #shows FFT from periodogram function for both aspirated length as pressure
REF = True             #shows isolated oscillation with most comparable 
Emod = True            #shows elastic modulus of single measurement for both FFT or Lock-In determination
Eav = True             #shows average elastic modulus of multiple measurements
A_P_av = True          #shows average amplitude of aspirated length and phase difference between pressure ans aspirated length signal

#VARIABLES THAT NEED TO BE GIVEN AS AN INPUT
amount_freqs = 7        #to determine how many frequencies have to be found in oscillation
Fs = 4000               #Sampling frequency
xN = 1                  #variable for how many times the data should be repeated by the tile function
N = 2                   #how many cycles for the lowest frequency
order = 2               #order of lowpass filter in lock-in

#define phase and amplitude
amplitude = [50]*amount_freqs
phase = [0]*amount_freqs

#data for generation fake data
frequencies = [0.05, 0.1, 0.25, 0.5, 1, 2, 6]
amplitude2 = [30]*amount_freqs
phase2 = [1/4*math.pi]*amount_freqs

#radius pipette and sample
rad_pip = 43e-6
rad_sam = 156e-6

#Counting amount of data files
amount = 0

#lists for finding the average elastic modulus, amplitude (aspirated length) and phase lag
E_storage_tot = [0] * amount_freqs          #save Estorage for single measurement determined with Lock-in
E_loss_tot = [0] * amount_freqs             #save Eloss for single measurement determined with Lock-in
Tand_tot = [0] * amount_freqs               #save tangend for single measurement determined with Lock-in
E_storage_tot_FFT = [0] * amount_freqs      #save Estorage for single measurement determined with FFT
E_loss_tot_FFT = [0] * amount_freqs         #save Eloss for single measurement determined with FFT
Tand_tot_FFT = [0] * amount_freqs           #save tangend for single measurement determined with FFT
AMP_tot = [0] * amount_freqs                #save amplitude for single measurement determined with Lock-in
PHA_tot = [0] * amount_freqs                #save phase for single measurement determined with Lock-in
AMP_tot_FFT = [0] * amount_freqs            #save amplitude for single measurement determined with FFT
PHA_tot_FFT = [0] * amount_freqs            #save phase for single measurement determined with FFT

lists_loss = {}                             #list of empty lists of E_loss per frequency from multiple measurements determined with Lock-in
lists_storage = {}                          #list of empty lists of E_storage per frequency from multiple measurements determined with Lock-in
lists_tand = {}                             #list of empty lists of tangend per frequency from multiple measurements determined with Lock-in
lists_loss_FFT = {}                         #list of empty lists of E_loss per frequency from multiple measurements determined with FFT
lists_storage_FFT = {}                      #list of empty lists of E_storage per frequency from multiple measurements determined with FFT
lists_tand_FFT = {}                         #list of empty lists of tangend per frequency from multiple measurements determined with FFT
lists_AMP = {}                              #list of empty lists of amplitude per frequency from multiple measurements determined with Lock-in
lists_PHA = {}                              #list of empty lists of phase per frequency from multiple measurements determined with Lock-in
lists_AMP_FFT ={}                           #list of empty lists of amplitude per frequency from multiple measurements determined with FFT
lists_PHA_FFT  = {}                         #list of empty lists of phase per frequency from multiple measurements determined with FFT

#create list inside the lists above for each frequency
for i in range (0,amount_freqs):
    lists_loss["f_loss"+str(i+1)] = []
    lists_storage["f_storage"+str(i+1)] = []
    lists_tand["f_tand"+str(i+1)] = []
    lists_loss_FFT["f_loss_FFT"+str(i+1)] = []
    lists_storage_FFT["f_storage_FFT"+str(i+1)] = []
    lists_tand_FFT["f_tand_FFT"+str(i+1)] = []
    lists_AMP["AMP"+str(i+1)] = []
    lists_PHA["PHA"+str(i+1)] = []
    lists_AMP_FFT["AMP_FFT"+str(i+1)] = []
    lists_PHA_FFT["PHA_FFT"+str(i+1)] = []

#loop over deltasense data that are in "Data_Eline", they should also be present in the file from which the code is running. In my case that is Python_code
for filename in os.listdir(str("Data_Eline")) :
    if filename.endswith (".tdms") :
        tdms_file = TdmsFile.read(filename)
        group = tdms_file['Demodulated data']
        p_channel = group["pressure"]
        l_channel = group["aspirated length"] 
        
        #filter data with lowpass, cutoff = 10 Hz
        p_channel = def_butter_lowpass.butter_lowpass_filter(p_channel, 10, Fs)
        l_channel = def_butter_lowpass.butter_lowpass_filter(l_channel, 10, Fs)

        #Define oscillation, preload, wait and unload are all filtered out.
        time, pressure, aspirated_length = def_period.define_period(p_channel, xN, l_channel, perio, Fs)           
        
        #define amount of operations
        operations = len(pressure)
        
        #generation of fake data to see how program works. 
        #time, pressure = def_reference.reference_data(frequencies, amplitude, phase, operations)
        #Stime, aspirated_length = def_reference.reference_data(frequencies, amplitude2, phase2, operations)

        #Define displayed frequencies by finding the maxima of the FFT. FFT defined by using periodogram function
        f, maxima_f, max_freqs, f, Pxx_den = def_FFT.maxima(pressure, amount_freqs, xN, N, Fs)
        max_frequency = max(max_freqs)
        min_frequency = min(max_freqs)
        
        f_L, maxima_f_L, max_freqs_L, f_L, Pxx_den_L = def_FFT.maxima(aspirated_length, amount_freqs, xN, N, Fs)
        max_frequency = max(max_freqs)
        min_frequency = min(max_freqs)

        #Finding amplitude and phase by using FFT. FFT defined by np.fft.fft function.
        fft_data_P, freq_P, magnitude2_P, phase2_P = def_FFT2.ampli_phase_FFT(pressure, max_freqs, time, Fs)
        fft_data_L, freq_L, magnitude2_L, phase2_L = def_FFT2.ampli_phase_FFT(aspirated_length, max_freqs, time, Fs)
        

        #lists of magintude, phase, phase difference, elastic storage and loss and the difference between both
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
        Tand_FFT = []

        #counter to call the right frequency from the list.
        x = 0

        #Filter out single frequencies
        for frequency in max_freqs:
            
            #find the amplitude and phase of the PRESSURE signal by lock-in principle
            R_freq_P, θ_freq_P = def_lock_in.filter_freq(pressure, frequency, max_frequency, order, Fs)

            #save amplitude and phase in lists
            R_P.append(R_freq_P)
            θ_P.append(θ_freq_P)

            #find the amplitude and phase of the ASPIRATED LENGTH signal by lock-in principle
            R_freq_L, θ_freq_L = def_lock_in.filter_freq(aspirated_length, frequency, max_frequency, order, Fs)
            R_L_freq = R_freq_L / (1* math.exp(9))  #multiply by e^-9 to have the amplitude in meters

            #save amplitude and phase in lists
            R_L.append(R_L_freq)
            θ_L.append(θ_freq_L)

            #calculate difference between phases and append to list
            θ = (θ_freq_P - θ_freq_L)
            θ_T.append(θ)

            #find elastic storage and loss of signal with usage of lock-in for A and θ determination
            E_S, E_L, T = def_E_modulus.elas_modulus(rad_pip, rad_sam, R_freq_P, R_L_freq, θ)

            #append lists with variables obtained from function above
            E_storage.append(E_S)
            E_loss.append(E_L)
            Tand.append(T)   

            #call amplitude and phase per frequency. x is counting which frequency should be called
            magni_FFT_P = magnitude2_P[x]
            magni_FFT_L = magnitude2_L[x] / (1* math.exp(9))
            phase_FFT_P = phase2_P[x]
            phase_FFT_L = phase2_L[x]

            #calculate difference between phases and append to list
            θ_FFT = (phase_FFT_P - phase_FFT_L)
            θ_T_FFT.append(θ_FFT)
            
            #find elastic storage and loss of signal with usage of FFT for A and θ determination
            E_S_FFT, E_L_FFT, T_FFT = def_E_modulus.elas_modulus(rad_pip, rad_sam, magni_FFT_P, magni_FFT_L, θ_FFT)  #θ_FFT

            E_storage_FFT.append(E_S_FFT)
            E_loss_FFT.append(E_L_FFT)
            Tand_FFT.append(T_FFT) 

            #counter for frequency calling
            x = x + 1

        #plot new amplitude and phase as a reference graph. A and θ determined from Lock-In
        ref_P_time, ref_P = def_reference.reference_data(max_freqs, R_P, θ_P, operations, Fs)
        ref_L_time, ref_L = def_reference.reference_data(max_freqs, R_L, θ_L, operations, Fs)

        #plot new amplitude and phase as a reference graph. A and θ determined from FFT
        ref_FFT_time_P, ref_FFT_P = def_reference.reference_data(max_freqs, magnitude2_P, phase2_P, operations, Fs)
        ref_FFT_time_L, ref_FFT_L = def_reference.reference_data(max_freqs, magnitude2_L, phase2_L, operations, Fs)

        #find coherence of lock-in  and FFT reference compared to the signal
        co_P = signal.coherence(ref_P, pressure, Fs, window = "hamming", nperseg = int(3*1/min_frequency))
        co_L = signal.coherence(ref_L, aspirated_length, Fs, window = "hamming", nperseg = int(3*1/min_frequency))
        co_P_FFT = signal.coherence(ref_FFT_P, pressure, Fs, window = "hamming", nperseg = int(3*1/min_frequency))
        co_L_FFT = signal.coherence(ref_FFT_L, aspirated_length, Fs, window = "hamming", nperseg = int(3*1/min_frequency))

        print(f"\nCoherence: \nPressure lock-in: {co_P}\nAspirated length lock-in: {co_L} \nPressure FFT: {co_P_FFT} \nAspirated length FFT: {co_L_FFT}\n")

        #print usefull information of a single measurement
        print(f"\nFrequencies pressure: {max_freqs}\n\nFrequencies aspirated length: {max_freqs_L}  \n\nMaginitude Pressure: {R_P} \nPhase Pressure: {θ_P} \n\nMaginitude Aspirated Lenght: {R_L} \nPhase Aspirated Length: {θ_L} \n\nPhase difference: {θ_T}\n \nE': {E_storage}\nE'': {E_loss} \nE''/E' = {Tand} \n")
        print(f"\nMagnitude FFT Pressure: {magnitude2_P} \nPhase FFT pressure: {phase2_P}\n\nMagnitude FFT Aspirated Length: {magnitude2_L} \nPhase FFT Aspirated Length: {phase2_L}\n\nPhase difference FFT: {θ_T_FFT} \n\nE' FFT: {E_storage_FFT} \nE'' FFT: {E_loss_FFT}  \nE''/E' FFT: {Tand_FFT} \n")

        #plot whole measurement for both aspirated length as pressure. Shows plot when real is True
        if real == True:

            fig4, ax1 = plt.subplots()
            ax1.plot(p_channel, label = "Pressure")
            ax1.set_xlabel("operations")
            ax1.set_ylabel("Pressure [Pa]")
            plt.legend(loc="upper left")
            ax2 = ax1.twinx()
            ax2.plot(l_channel, label = "Aspirated length", color = "#A60628")
            ax2.set_ylabel("Displacement [nm]")
            plt.legend(loc="upper right")

            fig5, ax1 = plt.subplots()
            ax1.plot(l_channel, -p_channel)
            ax1.set_xlabel("Displacement [nm]")
            ax1.set_ylabel("Pressure [Pa]")

        #plot of FFT from peiodogram function for both aspirated length as pressure. Shows plot when FFT is True
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

        #plot of oscillations including the reconstructed oscillations for both FFT and lock-in determination. Shows plot when REF is True
        if REF == True:
            plt.figure(8)
            plt.plot(time, pressure, label = "Signal")
            plt.plot(ref_P_time, ref_P, label = "Reference")
            plt.title("Pressure signal determined from Lock-In")
            plt.xlabel("Time [s]")
            plt.ylabel("Pressure [Pa]")
            plt.legend(loc="upper left")

            plt.figure(9)
            plt.plot(time, aspirated_length, label = "Signal")
            plt.plot(ref_L_time, ref_L, label = "Reference")
            plt.title("Aspirated length signal determined from Lock-In")
            plt.xlabel("Time [s]")
            plt.ylabel("δOPL [nm]")
            plt.legend(loc="upper left")

            plt.figure(10)
            plt.plot(time, pressure, label = "Signal")
            plt.plot(ref_FFT_time_P, ref_FFT_P, label = "Reference")
            plt.title("Pressure signal determined from FFT")
            plt.xlabel("Time [s]")
            plt.ylabel("Pressure [Pa]")
            plt.legend(loc="upper left")

            plt.figure(11)
            plt.plot(time, aspirated_length, label = "Signal")
            plt.plot(ref_FFT_time_L, ref_FFT_L, label = "Reference")
            plt.title("Aspirated length signal determined from FFT")
            plt.xlabel("Time [s]")
            plt.ylabel("δOPL [nm]")
            plt.legend(loc="upper left")

            fig12, ax1 = plt.subplots()
            ax1.plot(aspirated_length, pressure, label = "Signal")
            ax1.set_xlabel("Displacement [nm]")
            ax1.set_ylabel("Pressure [Pa]")
            plt.legend(loc="upper left")
            ax2 = ax1.twinx()
            ax2.plot(ref_FFT_L, ref_FFT_P, label = "Reference FFT", color = "#A60628")
            ax2.plot(ref_L, ref_P, label = "Reference Lock-In", color = "#7A68A6")
            ax2.set_ylabel("Pressure [Pa]")
            plt.legend(loc="upper right")
        
        #plot of elastic modulus of a single measurement     
        if Emod == True:
            fig13, ax1 = plt.subplots()
            ax1.loglog(max_freqs, E_storage, label = "E'")
            ax1.scatter(max_freqs, E_storage)
            ax1.loglog(max_freqs, E_loss, label = "E''")
            ax1.scatter(max_freqs, E_loss)
            ax1.set_xlabel("Frequency [Hz]")
            ax1.set_ylabel("E', E'' [kPa]")
            plt.legend(loc="upper left")
            
            ax2 = ax1.twinx()
            ax2.plot(max_freqs, Tand, label = "tan δ", color = "#7A68A6")
            ax2.scatter(max_freqs, Tand, color = "#7A68A6")
            ax2.set_ylabel("tan δ")
            plt.legend(loc="upper right")
            plt.title("Elastic moduli with Lock-In determination")
           
            fig14, ax1 = plt.subplots()
            ax1.loglog(max_freqs, E_storage_FFT, label = "E' ")
            ax1.scatter(max_freqs, E_storage_FFT)
            ax1.loglog(max_freqs, E_loss_FFT, label = "E''")
            ax1.scatter(max_freqs, E_loss_FFT)
            ax1.set_xlabel("Frequency [Hz]")
            ax1.set_ylabel("E', E'' [kPa]")
            plt.legend(loc="upper left")
            
            ax2 = ax1.twinx()
            ax2.plot(max_freqs, Tand_FFT, label = "tan δ", color = "#7A68A6")
            ax2.scatter(max_freqs,Tand_FFT, color = "#7A68A6")
            ax2.set_ylabel("tan δ")
            plt.legend(loc="upper right")
            plt.title(f"{filename}\nElastic moduli with FFT determination")

        #show plots of single measurement
        plt.show()

        #fill lists with E loss and storage and tangend
        for i in range(0, amount_freqs):
            lists_loss["f_loss"+str(i+1)].append(E_loss[i])
            lists_storage["f_storage"+str(i+1)].append(E_storage[i])
            lists_tand["f_tand"+str(i+1)].append(Tand[i])

            lists_loss_FFT["f_loss_FFT"+str(i+1)].append(E_loss_FFT[i])
            lists_storage_FFT["f_storage_FFT"+str(i+1)].append(E_storage_FFT[i])
            lists_tand_FFT["f_tand_FFT"+str(i+1)].append(Tand_FFT[i])

            lists_AMP["AMP"+str(i+1)].append(R_L[i])
            lists_PHA["PHA"+str(i+1)].append(θ_T[i])
            lists_AMP_FFT["AMP_FFT"+str(i+1)].append(magnitude2_L[i])
            lists_PHA_FFT["PHA_FFT"+str(i+1)].append(θ_T_FFT[i])
    
        #Sum all elastic moduli to eventually calculate the mean
        E_storage_tot = np.add(E_storage_tot, E_storage)
        E_loss_tot = np.add(E_loss_tot, E_loss)
        Tand_tot = np.add(Tand_tot, Tand)

        E_storage_tot_FFT = np.add(E_storage_tot_FFT, E_storage_FFT)
        E_loss_tot_FFT = np.add(E_loss_tot_FFT, E_loss_FFT)
        Tand_tot_FFT = np.add(Tand_tot_FFT, Tand_FFT)

        AMP_tot = np.add(AMP_tot, R_L)
        PHA_tot = np.add(PHA_tot, θ_T)
        AMP_tot_FFT = np.add(AMP_tot_FFT, magnitude2_L)
        PHA_tot_FFT = np.add(PHA_tot_FFT, θ_T_FFT)

        #count the amount of files
        amount = amount + 1

#calculate mean values    
E_storage_mean = E_storage_tot / amount
E_loss_mean = E_loss_tot / amount
Tand_mean = Tand_tot / amount

E_storage_mean_FFT = E_storage_tot_FFT / amount
E_loss_mean_FFT = E_loss_tot_FFT / amount
Tand_mean_FFT = Tand_tot_FFT / amount

AMP_mean = AMP_tot / amount
PHA_mean = PHA_tot / amount
AMP_mean_FFT = AMP_tot_FFT / amount
PHA_mean_FFT = PHA_tot_FFT / amount

#make lists for standard deviation
Eloss_std_dev = []
Estorage_std_dev = []
Tand_std_dev = []
Eloss_std_dev_FFT = []
Estorage_std_dev_FFT = []
Tand_std_dev_FFT = []
AMP_std_dev = []
PHA_std_dev = []
AMP_std_dev_FFT = []
PHA_std_dev_FFT = []

#Loop over frequencies to calculate standard deviation
for i in range(0, amount_freqs):
    #call lists with values per frequency
    std_dev_loss = lists_loss["f_loss"+str(i + 1)]
    std_dev_storage = lists_storage["f_storage"+str(i + 1)]
    std_dev_tand = lists_tand["f_tand"+str(i + 1)]
    std_dev_loss_FFT = lists_loss_FFT["f_loss_FFT"+str(i + 1)]
    std_dev_storage_FFT = lists_storage_FFT["f_storage_FFT"+str(i + 1)]
    std_dev_tand_FFT = lists_tand_FFT["f_tand_FFT"+str(i + 1)]
    std_dev_AMP = lists_AMP["AMP"+str(i+1)]
    std_dev_PHA = lists_PHA["PHA"+str(i+1)]
    std_dev_AMP_FFT = lists_AMP_FFT["AMP_FFT"+str(i+1)]
    std_dev_PHA_FFT = lists_PHA_FFT["PHA_FFT"+str(i+1)]

    #make list of absolute differences compared to the mean of 1 frequency for multiple measurements
    lst_std_dev_loss = []
    lst_std_dev_storage = []
    lst_std_dev_tand = []
    lst_std_dev_loss_FFT = []
    lst_std_dev_storage_FFT = []
    lst_std_dev_tand_FFT = []
    lst_std_dev_AMP = []
    lst_std_dev_PHA = []
    lst_std_dev_AMP_FFT = []
    lst_std_dev_PHA_FFT = []


    #loop over amount of measurements and determine absolute difference of measurements of certain frequency compared to the mean
    for x in range(0, amount):
        std_dev_L = (std_dev_loss[x] - E_loss_mean[i]) ** 2
        std_dev_S = (std_dev_storage[x] - E_storage_mean[i]) ** 2
        std_dev_T = (std_dev_tand[x] - Tand_mean[i]) ** 2
        std_dev_L_FFT = (std_dev_loss_FFT[x] - E_loss_mean_FFT[i]) ** 2
        std_dev_S_FFT = (std_dev_storage_FFT[x] - E_storage_mean_FFT[i]) ** 2
        std_dev_T_FFT = (std_dev_tand_FFT[x] - Tand_mean_FFT[i]) ** 2
        std_dev_A = (std_dev_AMP[x] - AMP_mean[i]) ** 2
        std_dev_P = (std_dev_PHA[x] - PHA_mean[i]) ** 2
        std_dev_A_FFT = (std_dev_AMP_FFT[x] - AMP_mean_FFT[i]) ** 2
        std_dev_P_FFT = (std_dev_PHA_FFT[x] - PHA_mean_FFT[i]) ** 2

        #add absolute difference of each measurement to lists
        lst_std_dev_loss.append(std_dev_L)
        lst_std_dev_storage.append(std_dev_S)
        lst_std_dev_tand.append(std_dev_T)
        lst_std_dev_loss_FFT.append(std_dev_L_FFT)
        lst_std_dev_storage_FFT.append(std_dev_S_FFT)
        lst_std_dev_tand_FFT.append(std_dev_T_FFT)
        lst_std_dev_AMP.append(std_dev_A)
        lst_std_dev_PHA.append(std_dev_P)
        lst_std_dev_AMP_FFT.append(std_dev_A_FFT)
        lst_std_dev_PHA_FFT.append(std_dev_P_FFT)


    #calculate final standard deviation of this frequency i
    std_dev_L_fin = math.sqrt(sum(lst_std_dev_loss) / amount)
    std_dev_S_fin = math.sqrt(sum(lst_std_dev_storage) / amount)
    std_dev_T_fin = math.sqrt(sum(lst_std_dev_tand) / amount)
    std_dev_L_fin_FFT = math.sqrt(sum(lst_std_dev_loss_FFT) / amount)
    std_dev_S_fin_FFT = math.sqrt(sum(lst_std_dev_storage_FFT) / amount)
    std_dev_T_fin_FFT = math.sqrt(sum(lst_std_dev_tand_FFT) / amount)
    std_dev_A_fin = math.sqrt(sum(lst_std_dev_AMP) / amount)
    std_dev_P_fin = math.sqrt(sum(lst_std_dev_PHA) / amount)
    std_dev_A_fin_FFT = math.sqrt(sum(lst_std_dev_AMP_FFT) / amount)
    std_dev_P_fin_FFT = math.sqrt(sum(lst_std_dev_PHA_FFT) / amount)

    #add standard deviation of frequency to list with std_dev of all frequencies
    Eloss_std_dev.append(std_dev_L_fin)
    Estorage_std_dev.append(std_dev_S_fin)
    Tand_std_dev.append(std_dev_T_fin)
    Eloss_std_dev_FFT.append(std_dev_L_fin_FFT)
    Estorage_std_dev_FFT.append(std_dev_S_fin_FFT)
    Tand_std_dev_FFT.append(std_dev_T_fin_FFT)
    AMP_std_dev.append(std_dev_A_fin)
    PHA_std_dev.append(std_dev_P_fin)
    AMP_std_dev_FFT.append(std_dev_A_fin_FFT)
    PHA_std_dev_FFT.append(std_dev_P_fin_FFT)

#print important variables
print(f"\nLOCK-IN \nEloss standard deviation: {Eloss_std_dev}  \nEstorage standard deviation: {Estorage_std_dev} \nTangend standard deviation: {Tand_std_dev}\n")
print(f"\nFFT \nEloss standard deviation: {Eloss_std_dev_FFT}  \nEstorage standard deviation: {Estorage_std_dev_FFT}\nTangend standard deviation: {Tand_std_dev_FFT}\n")

#plot average elastic moduli including standard deviation
if Eav == True:
    fig15, ax1 = plt.subplots()
    ax1.loglog(max_freqs, E_storage_mean, label = "E'", color = "#348ABD")
    ax1.scatter(max_freqs, E_storage_mean, color = "#348ABD")
    ax1.errorbar(max_freqs, E_storage_mean, yerr = Estorage_std_dev, color = "#348ABD")
    ax1.loglog(max_freqs, E_loss_mean, label = "E''", color = "#A60628")
    ax1.scatter(max_freqs, E_loss_mean, color = "#A60628")
    ax1.errorbar(max_freqs, E_loss_mean, yerr = Eloss_std_dev, color = "#A60628")
    ax1.set_xlabel("Frequency [Hz]")
    ax1.set_ylabel("E', E'' [kPa]")
    plt.legend(loc="upper left")
            
    ax2 = ax1.twinx()
    ax2.plot(max_freqs, Tand_mean, label = "tan δ", color = "#7A68A6")
    ax2.scatter(max_freqs, Tand_mean, color = "#7A68A6")
    ax2.errorbar(max_freqs, Tand_mean, yerr = Tand_std_dev, color = "#7A68A6")
    ax2.set_ylabel("tan δ")
    plt.legend(loc="upper right")
    plt.title(f"Average elastic modulus on {amount} samples determined with Lock-In")

    fig16, ax1 = plt.subplots()
    ax1.loglog(max_freqs, E_storage_mean_FFT, label = "E'", color = "#348ABD")
    ax1.scatter(max_freqs, E_storage_mean_FFT, color = "#348ABD")
    ax1.errorbar(max_freqs, E_storage_mean_FFT, yerr = Estorage_std_dev_FFT, color = "#348ABD")
    ax1.loglog(max_freqs, E_loss_mean_FFT, label = "E''", color = "#A60628")
    ax1.scatter(max_freqs, E_loss_mean_FFT, color = "#A60628")
    ax1.errorbar(max_freqs, E_loss_mean_FFT, yerr = Eloss_std_dev_FFT, color = "#A60628")
    ax1.set_xlabel("Frequency [Hz]")
    ax1.set_ylabel("E', E'' [kPa]")
    plt.legend(loc="upper left")
            
    ax2 = ax1.twinx()
    ax2.plot(max_freqs, Tand_mean_FFT, label = "tan δ", color = "#7A68A6")
    ax2.scatter(max_freqs, Tand_mean_FFT, color = "#7A68A6")
    ax2.errorbar(max_freqs, Tand_mean_FFT, yerr = Tand_std_dev_FFT, color = "#7A68A6")
    ax2.set_ylabel("tan δ")
    plt.legend(loc="upper right")
    plt.title(f"Average elastic modulus on {amount} samples determined with FFT")

    #plot for autocorrelation
    corr = signal.correlate(ref_L,aspirated_length, mode = "full")
    lags = signal.correlation_lags(len(ref_L), len(aspirated_length), mode = "full")
    corr /= np.max(corr)
    plt.figure(17)
    plt.plot(lags,corr)

if A_P_av == True:
    fig18, ax1 = plt.subplots()
    ax1.loglog(max_freqs, AMP_mean, label = "Aspirated length", color = "#348ABD")
    ax1.scatter(max_freqs, AMP_mean, color = "#348ABD")
    ax1.errorbar(max_freqs, AMP_mean, yerr = AMP_std_dev, color = "#348ABD")
    ax1.set_xlabel("Frequency [Hz]")
    ax1.set_ylabel("δOPL [nm]")
    plt.legend(loc="upper left")
            
    ax2 = ax1.twinx()
    ax2.plot(max_freqs, PHA_mean, label = "Phase lag", color = "#7A68A6")
    ax2.scatter(max_freqs, PHA_mean, color = "#7A68A6")
    ax2.errorbar(max_freqs, PHA_mean, yerr = PHA_std_dev, color = "#7A68A6")
    ax2.set_ylabel("Phase lag")
    plt.legend(loc="upper right")
    plt.title(f"Average amplitude and phase lag on {amount} samples determined with Lock-In")

    fig16, ax1 = plt.subplots()
    ax1.loglog(max_freqs, AMP_mean_FFT, label = "Amplitude", color = "#348ABD")
    ax1.scatter(max_freqs, AMP_mean_FFT, color = "#348ABD")
    ax1.errorbar(max_freqs, AMP_mean_FFT, yerr = AMP_std_dev_FFT, color = "#348ABD")
    ax1.set_xlabel("Frequency [Hz]")
    ax1.set_ylabel("δOPL [nm]")
    plt.legend(loc="upper left")
            
    ax2 = ax1.twinx()
    ax2.plot(max_freqs, PHA_mean_FFT, label = "Phase lag", color = "#7A68A6")
    ax2.scatter(max_freqs, PHA_mean_FFT, color = "#7A68A6")
    ax2.errorbar(max_freqs, PHA_mean_FFT, yerr = PHA_std_dev_FFT, color = "#7A68A6")
    ax2.set_ylabel("Phase lag")
    plt.legend(loc="upper right")
    plt.title(f"Average amplitude and phase lag on {amount} samples determined with FFT")

#show plot of average elastic moduli
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


