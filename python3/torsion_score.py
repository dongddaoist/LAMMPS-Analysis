from numpy import *
import math
import os.path
import numpy as np
import pandas as pd
from itertools import islice
from scipy.fft import fft, ifft,fftfreq

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import make_interp_spline,BSpline
import matplotlib.patches as m_patches


ifile="run2-shrink/torsion.dat"
ntort=58
head_label=['indx','angle']
iline=ntort+1
#list_of_rows=linspace(0,nframe,nframe+1)
df=pd.read_csv(ifile,skiprows=[0],keep_default_na=False,header=None)
nrows=df.shape[0]
nf=int(ceil(nrows/iline))

#print(df.shape)

rows_comm=np.array(linspace(0,nf-1,nf)*iline) #+1
r_int=rows_comm.astype(int)
box=df.iloc[r_int]
dataset=df.drop(df.index[r_int])
sd=dataset[0].str.split(expand=True)
sd.columns=head_label
sd[["indx","angle"]]=sd[["indx","angle"]].apply(pd.to_numeric)

#print(dataset.shape)
#print(box.shape)
#print(dataset.iloc[1:10].values)

indx_array=sd[["indx"]].values
angle_array=sd[["angle"]].values

index=indx_array.reshape((-1,ntort),order='C').astype(int)
angle=angle_array.reshape((-1,ntort),order='C')
#print(index.shape)

print(angle)
angle_ref=angle*0.0

#angle2=angle0*0.0
angle_ref[np.logical_and(index==5, angle<=90.0)]=60.0
angle_ref[np.logical_and(index==5, angle>90.0)]=120.0
angle_ref[index==7]=180.0
angle_ref[index==11]=180.0
#print(index)
#print(angle_ref)

print(angle.shape)

diff=abs(np.sin((angle-angle_ref)/180.0*np.pi))

f_score=np.sum(diff,axis=1)
print(f_score)
x=np.linspace(0.15,float(nf)*0.15,nf)*0.5
plt.plot(x,f_score,linewidth=2.0)
plt.show()
