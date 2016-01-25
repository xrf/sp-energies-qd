'''
R^-c extrapolation script
'''

from scipy.optimize import curve_fit
import numpy as np
import scipy.special as ss

def func(x,a,b,c):
    z=np.zeros(len(x))
    for i in range (len(x)):
        res = 0               
        for j in range (1, int(x[i])+1):
            res += j**(b)
        z[i] = c-res*a
    return z

'''
example:
@x shells
@y energies
@err error estimates
'''
#init arrays
x = np.array([5, 7, 9, 11, 13, 15, 17, 19, 21, 23]) #shells
y = np.array([15.016034, 14.990082, 14.978221, 14.971254, 14.967 , 14.963577, 14.961598 ,14.958911 , 14.958189, 14.957]) #energies
err = np.array([.001, .0005, .0005, .00047, .0009, .0006, .0006, .0008, .00055, .0007])

#fit curve
maxfev = 10000 #max iterations
gtol = 0.00001 #tolerance (see scipy manual)
popt, pcov = curve_fit(func, x, y, sigma=None, maxfev=10000, gtol=.00001)

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
