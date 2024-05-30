from pymatreader import read_mat
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spiceypy as sp
import datetime as dt
from tqdm import tqdm

mat = read_mat('../MATLAB/LHPOP/ORBdata.mat')
XJ2000 = mat["orb"]["XJ2000"]
XDF = pd.DataFrame(XJ2000,columns = ["X","Y","Z","VX","VY","VZ"])
sp.furnsh("lrorg_2021074_2021166_v01.bsp")
sp.furnsh("naif0012.tls")
start,stop = mat['orb']['epoch']['et']

### plot moon and traj
fig1 = plt.figure()
ax = plt.axes(projection='3d')

#plot moon
rM = 1738
u = np.linspace(0, 2 * np.pi, 50)
v = np.linspace(0, np.pi, 50)
x = rM * np.outer(np.cos(u), np.sin(v))
y = rM * np.outer(np.sin(u), np.sin(v))
z = rM * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z,cmap = "gray",zorder=1)
ax.set_aspect('equal')

#plot traj
l=None
ax.plot(XJ2000[:l,0],XJ2000[:l,1],XJ2000[:l,2],"r",zorder=10,label="traj all")
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')

# Plot initial position
state = sp.spkezr("LRO",start,"J2000","NONE","MOON")[0]
ax.plot(state[0],state[1],state[2],"b.",zorder=10,label="Initial pos")        #Starting point

# Plot final position
state = sp.spkezr("LRO",stop,"J2000","NONE","MOON")[0]
ax.plot(state[0],state[1],state[2],"c.",zorder=10,label="Final Pos True")        #Ending point True

# Plot final position Precise Prop
EndP = XJ2000[-1]
ax.plot(EndP[0],EndP[1],EndP[2],"y.",zorder=10,label="Final Pos All")        #Starting point All
ax.legend()
###

### Second graph
# acquire data
Tstep = mat['orb']['epoch']['span'][1] - mat['orb']['epoch']['span'][0]
Tspan = stop-start
Err = []
for _ in range(6):Err.append([])

for i in tqdm(range(round(Tspan/Tstep)+1)):
    Time = start + Tstep*i
    state = sp.spkezr("LRO",Time,"J2000","NONE","MOON")[0]
    for j in range(6):
        Err[j].append(state[j]-XJ2000[i,j])

# plotting the 6 errors graphs
fig2 = plt.figure(figsize=(18,4.8))
axes = fig2.subplots(2,3)
for i in range(6):
    axes[i//3,i%3].plot(np.arange(0,Tspan+Tstep/2,Tstep),Err[i],
                        label="RMSE = "+'{:.2e}'.format(np.sqrt((np.array(Err[i])**2).sum())/len(Err[i]))+[" km"," km/s"][i//3])
    axes[i//3,i%3].set_xlabel('s')
    axes[i//3,i%3].set_ylabel(f'error in {["X (km)","Y (km)","Z (km)","VX (km/s)","VY (km/s)","VZ (km/s)"][i]}')
    axes[i//3,i%3].plot([0,Tspan],[0,0])
    axes[i//3,i%3].legend()

# close everything and plot
sp.kclear()
plt.show()

# sp.dafec(sp.spklef("lrorg_2021074_2021166_v01.bsp"),1000)