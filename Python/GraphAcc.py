import matplotlib.pyplot as plt
from pymatreader import read_mat
import numpy as np

a=1

if a: file="acc"
else: file="alt"
with open(f"../MATLAB/LHPOP/output/{file}.txt") as file:
    Lines = file.readlines()

Acc = dict()
nv = 9
n=len(Lines)//(nv+1)
for i in range(n):
    alt=float(Lines[(nv+1)*i].strip().split()[1])
    Acc[alt]=[]
    for j in range(nv):
        data = Lines[(nv+1)*i+j+1].strip().split()
        Acc[alt].append(float(data[1]))


### Plot
name = ['Sun', 'Jupiter', 'Moon', 'LGF', 'Earth', 'EGF',  'RelatCorr', 'SunRadPress', 'EarthAlb']
values = np.array(list(Acc.values()))
fig,ax = plt.subplots()
for i in [2,4,3,0,7,6,5,1,8]:
    t=np.array(list(Acc.keys()))-list(Acc.keys())[0]
    ax.plot(t/3600,np.log10(values[:,i]),label=name[i])
if not a:
    ax.set_xlabel('Distance from moon\'s center (km)')
    ax.set_ylabel('Acceleration (log10) (km/s^2)')
    ax.set_title("Part of each perturbations in the total acceleration with respect to the altitude")
else:
    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('Acceleration (log10) (km/s^2)')
    ax.set_title("Part of each perturbations in the total acceleration with respect to the time")

ax.legend()