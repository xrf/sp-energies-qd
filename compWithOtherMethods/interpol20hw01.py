from pylab import *
import sys
import os
import numpy as np
from scipy.optimize import curve_fit

def func(x, a, b, c):
     return a+ b*(x**(c))
#     return a+ b*exp(-x*c)

dataSRG = loadtxt("20part01/20part.dat")
x = dataSRG[:,0]
y=dataSRG[:,2]

#popt, pcov = curve_fit(func, x[-4:-1], y[-4:-1], sigma=None, maxfev = 200000, gtol=0.00001)

popt, pcov = curve_fit(func, x,y, sigma=None, maxfev = 200000, gtol=0.00001)

a = popt[0]
b = popt[1]
c = popt[2]
print a
print b
print c
