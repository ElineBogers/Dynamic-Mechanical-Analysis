import matplotlib.pyplot as plt
from nptdms import TdmsFile
import csv
import math

def variance_cal(data_file):

    #number of operations
    n = len(data_file)

    # Mean of the data
    mean = sum(data_file) / n

    # Square deviations
    deviations = [(x - mean) ** 2 for x in data_file]
    
    # standard deviation
    standard_deviation = 2 * math.sqrt(sum(deviations) / n)
    
    return standard_deviation, mean



