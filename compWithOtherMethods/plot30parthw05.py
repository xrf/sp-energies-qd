#!/usr/bin/python
from numpy import *
from pylab import *
from matplotlib import rc, rcParams
import matplotlib.units as units
import matplotlib.ticker as ticker
import sys
import os


rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['New fits of nuclear forces']})

dataSRG = loadtxt("30part05/30part.dat")
dataCCSD = loadtxt("30part05/CCSD.dat")

l = len(dataSRG[:,0])
data7 = zeros(l)
data7[:] = 187.0426

axis([7,20,185, 202])
xlabel(r'$R$',fontsize=20)
ylabel(r'$E_0\, [E_h^*]$',fontsize=20)
plot(dataSRG[:,0], dataSRG[:,2],'b-*', dataSRG[:,0], dataSRG[:,1],'r:.', dataCCSD[:,0], dataCCSD[:,1], 'm:v',dataSRG[:,0], data7[:], 'k--',markersize=7)

#title(r'$N=30,\;\omega=0.5$', fontsize=20, horizontalalignment='center')
legend(('SRG','HF','CCSD','DMC'),
           'upper right', shadow=False, fancybox=False,prop={"size":18})
legend(loc='upper right')

xticks( [7,8,9,10,11,12,13,14,15,16,17,18,19,20])
# Save the figure in a separate file
savefig('30parthw05.pdf', format='pdf')
# Draw the plot to screen
show()
    
