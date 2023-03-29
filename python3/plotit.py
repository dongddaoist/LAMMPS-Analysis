import numpy  as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
nskip=0

data = np.loadtxt('acf-tot.dat',skiprows=nskip)
levels=np.linspace(0.0,0.5,21)

x = data[:,0]
y = data[:, 1]
lenx=len(x)
leny=len(y)

plt.plot(x, y,'b-')
plt.xlabel('$t\ $ [ps]',fontsize=14)
plt.ylabel('$ACF$',fontsize=14)
plt.xscale("log")
plt.yscale("log")
#plt.legend()
plt.show()

