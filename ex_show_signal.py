import matplotlib.pyplot as plt
from nptdms import TdmsFile
import def_period
import numpy as np
import os


for filename in os.listdir(str("Data_Eline")) :
    if filename.endswith (".tdms") :
        tdms_file = TdmsFile.read(filename)
        group = tdms_file['Demodulated data']
        p_channel = group["pressure"]
        p_channel_data = p_channel[:]
        l_channel = group['pressure']
        l_channel_data = l_channel[:]

        fig, ax1 = plt.subplots()
        ax1.plot(p_channel_data, label = "Pressure", color = 'royalblue')
        ax1.set_xlabel("operations")
        ax1.set_ylabel("Pressure [Pa]")

        ax2 = ax1.twinx()
        ax2.plot(l_channel_data, label = "Displacement", color = "firebrick")
        ax2.set_ylabel("Displacement [nm]")
        
        plt.legend(loc="upper left")

        plt.show()

 


#time, oscillation = defineperiod.define_period(p_channel_data))

#a_file = open("test.txt", "w")
#for row in oscillation:
    #np.savetext(a_file, row)
#a_file.close()

#plt.plot(time, oscillation)


