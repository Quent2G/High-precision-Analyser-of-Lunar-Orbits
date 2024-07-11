import matplotlib.pyplot as plt
from pymatreader import read_mat
import numpy as np

mat = read_mat('../MATLAB/LHPOP/output/ORBdata.mat')
T = mat["orb"]["t"]
with open("../MATLAB/LHPOP/output/acc.txt") as file:
    Lines = file.readlines()

Acc = { "SunMass":[],
        "JupMass":[],
        "LunGrav":[],
        "EarthGrav":[],
        "RelatCorr":[],
        "SunRadPress":[],
        "EarthAlb":[]}

n=len(Lines)//8
for i in range(n):
    t=float(Lines[8*i].strip().split()[1])
    for j in range(7):
        data = Lines[8*i+j+1].strip().split()
        Acc[list(Acc.keys())[j]].append((t,float(data[1])))

### separately
fig,ax = plt.subplots(2,4)
for i in range(7):
    data = np.array(Acc[list(Acc.keys())[i]])
    ax[i//4,i%4].plot((data[:,0]-T[0])/3600,data[:,1])
    ax[i//4,i%4].set_title(list(Acc.keys())[i])
    ax[i//4,i%4].set_xlabel('hours')
    ax[i//4,i%4].set_ylabel('acc')

### Altogether
fig,ax = plt.subplots()
for i in range(7):
    data = np.array(Acc[list(Acc.keys())[i]])
    ax.plot((data[:,0]-T[0])/3600,np.log10(data[:,1]),label=list(Acc.keys())[i])
ax.set_xlabel('hours')
ax.set_ylabel('acc (10^)')
ax.legend()
