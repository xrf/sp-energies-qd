from pylab import *
import sys
import os
import scipy.special as ss
import numpy as np
from scipy.optimize import curve_fit

def func(x,a,b,c):
    z=np.zeros(len(x))
    for i in range (len(x)):
        res = 0 
        for j in range (1, int(x[i])+1):
            res += j**(b)
        z[i] = c-res*a
    return z

dataSRG = loadtxt("6part10/FCI6part.dat")
x = dataSRG[:,0]
y = dataSRG[:,1]



popt, pcov = curve_fit(func, x[-4:-1], y[-4:-1], sigma=None, maxfev=200000, gtol=.00001)
a = popt[2]
b = popt[0]
c = -popt[1]
extrapolated_energy = a-b*ss.zeta(c, 1)
print extrapolated_energy
print a
print b
print c
