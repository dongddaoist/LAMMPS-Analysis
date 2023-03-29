from numpy import *
from scipy.optimize import curve_fit
from scipy import integrate
from scipy.integrate import quad
from function import funct
import numpy as np
import matplotlib.pyplot as plt


nskip=0
dstart=2
dend=200


data = np.loadtxt('acf-tot.dat',skiprows=nskip)
x0=data[:,0]
y0=data[:,1]
x1=data[dstart:dend,0]
y1=data[dstart:dend,1]
popt,pcov=curve_fit(funct,x1,y1,bounds=([0.7,0.2,1],[1.0,1.0,1000000]))
print(popt)

a=popt[0]
b=popt[1]
c=popt[2]
def  intef(x,a,b,c):
    return a * np.exp(-(x/c)**b)

res,err = quad(intef,0,100000,args=(a,b,c))
print(res)
xdata=linspace(10,10000,50000)
plt.plot(x0,y0,label='Original-data')
plt.plot(x0, funct(x0, *popt),label='Fitted')
plt.xlabel('$t\ $ [ps]',fontsize=14)
plt.ylabel('$ACF$',fontsize=14)
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.show()
