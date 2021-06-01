import math
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile
import def_std_dev
import def_preload
import def_start
import def_unload
import def_stop
import def_creep_equ
import pandas as pd

#within this piece of code the period of the oscillation of the DMA is determined. The data consists out of:
    #1 flat starting point (time it takes to start experiment, usually 1-3 seconds)
    #2 Preload : negative pressure applied in a linear manner (based on preload of 6-8 seconds, can also be more)
    #3 Wait : static period to give the sample some time to stabilize, should last at least 1 second, preferably more
    #4 oscillation : should at least last 1 second
    #5 wait : static period to give sample time to stabilize, should last at least 1 second
    #6 unload : positive pressure is applied linearly (based on unload of 6-8 seconds, can also be more)
    #7 flat stopping point (time it takes to stop recording)

#function that loops over data in file to find start
def define_period(pressure_data, N, length_data, perio, Fs) :

    #define standaard deviation of first part, take 0.75 seconds as interval. 
    sec_1 = int(0.75*Fs)
    data_std_dev = pressure_data[0:sec_1]
    std_dev_first, mean_first = def_std_dev.variance_cal(data_std_dev)

    #FILTER PRELOAD AND DETERMINE STANDARD DEVIATION
    N_wait, N_wait_end = def_preload.define_wait(pressure_data, std_dev_first, Fs)

    #define data for standard deviation
    wait_data = pressure_data[N_wait:N_wait_end]

    #define standard deviation
    standard_deviation, mean = def_std_dev.variance_cal(wait_data)

    #define starting point
    N_start = def_start.start_osci(N_wait, mean, std_dev_first, pressure_data, Fs)

    #define unload by reversing data
    N_wait_reversed, N_normal, N_wait_reversed_start = def_unload.define_unload(pressure_data, std_dev_first, Fs)

    #define data for reversed standard deviation
    rev_wait_data = pressure_data[N_wait_reversed_start:N_wait_reversed]

    #define standard deviation
    rev_std_dev, rev_mean = def_std_dev.variance_cal(rev_wait_data)
    
    #define stopping point
    N_stop = def_stop.stop_osci(N_wait_reversed, N_normal, rev_std_dev, rev_mean, pressure_data, Fs)

    #find mean of oscillation
    P_data_mean = pressure_data[N_start:N_stop]
    P_fin_std_dev, P_fin_mean = def_std_dev.variance_cal(P_data_mean)

    #define data oscillation
    data_pressure = ((pressure_data[N_start:N_stop]) - P_fin_mean) *76.59  
    len_oscillation = len(data_pressure)  
    data_pressure = np.tile(data_pressure, N)

    #amount of operations
    N_difference = N_stop - N_start

    #define starting time
    time = 0

    #list for timesteps and pressure steps
    time_definition = [] * N_difference 
    define_pressure = [] * N_difference 
    define_length = [] * N_difference

        
    #loop over time to define time steps
    for pressure in data_pressure:

        time = time + 1/Fs

        define_pressure.append(pressure*-1)
        time_definition.append(time)

    #define parameters of function creep_data
    P0 = - pressure_data[N_wait]
    L0 = length_data[N_wait]

    #Define creep period and isolate data.
    creep_data = length_data[N_wait:N_start] - L0
    operations_creep = len(creep_data)

    #create time_domain of the same length as the creep data
    time_creep = []
    for x in range(0, operations_creep): 
        time = x * 1/Fs
        time_creep.append(time)

    #define period of wait + oscillation    
    tot_time = (len_oscillation + operations_creep) * 1/Fs

    #make fit for creep
    creep_data2 = def_creep_equ.remove_creep(time_creep, creep_data, P0, tot_time, Fs)
    subtract_data = creep_data2[operations_creep:]
    
    #create aspirated length oscillation data
    if all(length_data) != 0:

        #subtract creepfit from oscillation data
        L_data_raw = length_data[N_start:N_stop]
        L_data_mean = L_data_raw - subtract_data
        L_fin_std_dev, L_fin_mean = def_std_dev.variance_cal(L_data_mean)

        #define data length
        data_length = ((L_data_mean) - L_fin_mean) / 1.33 
        data_length = np.tile(data_length, N)
    
        for length in data_length:

            define_length.append(length)
    else:
        define_length = 0

    #define each period area
    area1 = pressure_data[0:750]
    op1 = []
    time1 = 0
    for x in range(0,750): op1.append(time1); time1 = time1 + 0.001
    area2 = pressure_data[4000:N_wait]
    op2 = []
    time2 = 0.001 * 4000
    for x in range(4000,N_wait): op2.append(time2); time2 = time2 + 0.001
    area3 = pressure_data[N_wait:N_start]
    op3 = []
    time3 = 0.001 * N_wait
    for x in range(N_wait, N_start): op3.append(time3); time3 = time3 + 0.001
    area4 = pressure_data[N_start:N_stop]
    op4 = []
    time4 = 0.001 * N_start
    for x in range(N_start, N_stop): op4.append(time4); time4 = time4 + 0.001
    area5 = pressure_data[N_stop:N_wait_reversed]
    op5 = []
    time5 = 0.001 * N_stop
    for x in range(N_stop, N_wait_reversed): op5.append(time5); time5 = time5 + 0.001
    
    #define time array for plot
    operations_tot = len(pressure_data)
    tot_time = 0.001 * operations_tot
    time_step_tot = np.linspace(0, tot_time, num = operations_tot)

    #plot whole measurement including the segmentation areas
    if perio == True:
        plt.figure(1)
        plt.plot(time_step_tot, pressure_data, label = "Signal")
        plt.plot(op1, area1, label = "Area 1")
        plt.plot(op2, area2, label = "Area 2")
        plt.plot(op3, area3, label = "Area 3")
        plt.plot(op4, area4, label = "Area 4")
        plt.plot(op5, area5, label = "Area 5")
        plt.xlabel("Time [s]")
        plt.ylabel("Pressure [Pa]")
        plt.legend(loc="upper right")

        #plot creep against creepfit
        plt.figure(2)
        plt.plot(creep_data2, label = "Creep reference")
        plt.plot(creep_data, label = "Creep Signal")
        plt.xlabel("Operations")
        plt.ylabel("δOPL [nm]")
        plt.legend(loc="upper right")
        plt.title("Creep fit")

        #plot corrected, uncorrected and creepfit data between oscillation interval
        plt.figure(3)
        plt.plot(L_data_raw, label = "Raw data")
        plt.plot(L_data_mean, label = "Corrected data")
        plt.plot(subtract_data, label = "Creep")
        plt.xlabel("Operations")
        plt.ylabel("δOPL [nm]")
        plt.title("Data correction")
        plt.legend(loc="upper right")

    return time_definition, define_pressure, define_length



            

