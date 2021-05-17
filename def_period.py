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

#function that loops over data in file to find start
def define_period(pressure_data, N, length_data) :

    #define standaard deviation of first part
    data_std_dev = pressure_data[0:750]
    std_dev_first, mean_first = def_std_dev.variance_cal(data_std_dev)
    #print(f" \nStd_begin: {std_dev_first}")

    #FILTER PRELOAD AND DETERMINE STANDARD DEVIATION
    N_wait, N_wait_end = def_preload.define_wait(pressure_data, std_dev_first)

    #define data for standard deviation
    wait_data = pressure_data[N_wait:N_wait_end]

    #define standard deviation
    standard_deviation, mean = def_std_dev.variance_cal(wait_data)
    #print(f"Std_creep: {standard_deviation}")

    #define starting point
    N_start = def_start.start_osci(N_wait, mean, std_dev_first, pressure_data)

    #define std dev end
    tot_operations = len(pressure_data)
    start_std_dev_stop = tot_operations - 1000
    data_std_dev2 = pressure_data[start_std_dev_stop:tot_operations]
    std_dev_last, mean_last = standard_deviation, mean = def_std_dev.variance_cal(data_std_dev2)
    #print(f"Std_end: {std_dev_last}")

    #define unload by reversing data
    N_wait_reversed, N_normal, N_wait_reversed_start = def_unload.define_unload(pressure_data, std_dev_first)

    #define data for reversed standard deviation
    rev_wait_data = pressure_data[N_wait_reversed_start:N_wait_reversed]

    #define standard deviation
    rev_std_dev, rev_mean = def_std_dev.variance_cal(rev_wait_data)
    #print(f"Std_creep_end: {rev_std_dev} \n")
    #define stopping point
    N_stop = def_stop.stop_osci(N_wait_reversed, N_normal, rev_std_dev, rev_mean, pressure_data)

    #print(f"N_stop - N_start: {N_stop - N_start}, N_start: {N_start}, N_stop: {N_stop}, N_wait: {N_wait}")

    #find mean of oscillation
    P_data_mean = pressure_data[N_start:N_stop]
    P_fin_std_dev, P_fin_mean = def_std_dev.variance_cal(P_data_mean)

    #define data oscillation
    data_pressure = ((pressure_data[N_start:N_stop]) - P_fin_mean) *76.59    
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

        time = time + 0.001

        define_pressure.append(pressure*-1)
        time_definition.append(time)



    if all(length_data) != 0:
        L_data_mean = length_data[N_start:N_stop]
        L_fin_std_dev, L_fin_mean = def_std_dev.variance_cal(L_data_mean)

        #define data length
        data_length = ((length_data[N_start:N_stop]) - L_fin_mean) / 1.33 
        data_length = np.tile(data_length, N)
    
        for length in data_length:

            define_length.append(length)
    else:
        define_length = 0

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

    operations_tot = len(pressure_data)
    tot_time = 0.001 * operations_tot
    time_step_tot = np.linspace(0, tot_time, num = operations_tot)
    
    plt.plot(time_step_tot, pressure_data, label = "Signal", color = "Black")
    plt.plot(op1, area1, label = "Area 1", color = "Red")
    plt.plot(op2, area2, label = "Area 2", color = "Yellow")
    plt.plot(op3, area3, label = "Area 3", color = "Green")
    plt.plot(op4, area4, label = "Area 4", color = "royalblue")
    plt.plot(op5, area5, label = "Area 5", color = "Orange")
    
    plt.xlabel("Time [s]")
    plt.ylabel("Pressure [Pa]")
    plt.legend(loc="upper right")
    plt.show()

    #define parameters of function creep_data
    P0 = pressure_data[N_wait]
    creep_data = length_data[N_wait:N_start]
    operations = len(creep_data)
    time_end = 0.001 * operations
    time_series = np.linspace(0, time_end, num = operations)
    params, params_covariance = def_creep_equ.remove_creep(time_series, creep_data, P0)

    return time_definition, define_pressure, define_length



            

