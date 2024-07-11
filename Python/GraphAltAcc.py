import matplotlib.pyplot as plt
from pymatreader import read_mat
import numpy as np

with open("../MATLAB/LHPOP/output/alt.txt") as file:
    Lines = file.readlines()

Acc = dict()

n=len(Lines)//7
for i in range(n):
    alt=float(Lines[7*i].strip().split()[1])
    Acc[alt]=[]
    for j in range(6):
        data = Lines[7*i+j+1].strip().split()
        Acc[alt].append(float(data[1]))


### Plot
name = ['SunMass', 'EarthMass', 'LunGrav', 'RelatCorr', 'SunRadPress', 'EarthAlb']
values = np.array(list(Acc.values()))
fig,ax = plt.subplots()
for i in range(6):
    data = np.array(Acc[list(Acc.keys())[i]])
    ax.plot(Acc.keys(),np.log10(values[:,i]),label=name[i])
ax.set_xlabel('Distance from moon\'s center (km)')
ax.set_ylabel('Acceleration (log10) (km/s^2)')
ax.legend()
ax.set_title("Part of each perturbations in the total acceleration with respect to the altitude")
