import matplotlib.pyplot as plt
from pymatreader import read_mat
import numpy as np

mat = read_mat('../MATLAB/LHPOP/output/ORBdata.mat')
T = mat["orb"]["t"]
with open("../MATLAB/LHPOP/output/acc.txt") as file:
    Lines = file.readlines()

Acc = { "SunMass":[],
        "EarthMass":[],
        "LunGrav":[],
        "RelatCorr":[],
        "SunRadPress":[],
        "EarthAlb":[]}

n=len(Lines)//6
for i in range(n):
    for j in range(6):
        data = Lines[6*i+j].strip().split()
        Acc[list(Acc.keys())[j]].append((float(data[1]),float(data[2])))

### separately
# fig,ax = plt.subplots(2,3)
# for i in range(6):
#     data = np.array(Acc[list(Acc.keys())[i]])
#     ax[i//3,i%3].plot((data[:,0]-T[0])/3600,data[:,1])
#     ax[i//3,i%3].set_title(list(Acc.keys())[i])
#     ax[i//3,i%3].set_xlabel('hours')
#     ax[i//3,i%3].set_ylabel('acc')

### Altogether
fig,ax = plt.subplots()
for i in range(6):
    data = np.array(Acc[list(Acc.keys())[i]])
    ax.plot((data[:,0]-T[0])/3600,np.log10(data[:,1]),label=list(Acc.keys())[i])
ax.set_xlabel('hours')
ax.set_ylabel('acc (10^)')
ax.legend()
