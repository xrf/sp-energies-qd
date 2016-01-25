import numpy as np
from scipy.optimize import curve_fit
import scipy.special as ss

def func(x,a,b,c):
    z=np.zeros(len(x))
    for i in range (len(x)):
        res = 0               
        for j in range (1, int(x[i])+1):
            res += j**(b)
        z[i] = c-res*a
    return z

#def func(x, a, b, c):
#     return a*(x**b) + c

x = np.array([10,12,14,16,18,20]) #shells
y = np.array([12.21380, 12.21819, 12.22080, 12.22152, 12.22196, 12.22229]) #energies

popt, pcov = curve_fit(func, x, y, sigma=None, maxfev = 10000, gtol=0.01)
#the slope is
slope = popt[1]
#the covariance matrix is
pcov

#the extrapolated energy is on closed form a-b*zeta(c)
#where ss.zeta is the Riemann-Zeta function
a = popt[2]
b = popt[0]
c = -popt[1]
extrapolated_energy = a-b*ss.zeta(c, 1)
print extrapolated_energy
print a
print b
print c
