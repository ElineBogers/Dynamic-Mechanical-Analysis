import math
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile
import os
from scipy.optimize import curve_fit

def fit_creepL(time, P0, K, n):
    return P0/K*(1 - np.exp(-(K/n*time)))


def remove_creep(time, data_aspirated, P0, tot_time):
    params, params_covariance = curve_fit(fit_creepL, time, data_aspirated, bounds = ([(P0-0.01), 0, 0], [(P0+0.01), np.inf, np.inf]))
    print(params)
    creep_data = []
    time_data = []

    time_mili = int(tot_time * 1000)
    for x in range(0, time_mili):

        if x < tot_time * 1000:
            x = x/1000
            creep = fit_creepL(x, P0, params[1], params[2])
            time_data.append(x)
            creep_data.append(creep)

    return creep_data
