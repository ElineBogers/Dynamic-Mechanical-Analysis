import math
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile
import os
from scipy.optimize import curve_fit

def fit_creepSL(time, P0, K, n):
    return P0/K*(1 - np.exp(-(K/n*time)))


def remove_creep(time, data_aspirated, P0):
    params, params_covariance = curve_fit(fit_creepSL, time, data_aspirated, bounds = ([(P0-0.01), 0, 0], [(P0+0.01), np.inf, np.inf]))
    return params, params_covariance
