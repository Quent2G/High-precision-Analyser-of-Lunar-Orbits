from pymatreader import read_mat
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spiceypy as sp
import datetime as dt
from tqdm import tqdm

mat = read_mat('../MATLAB/LHPOP/output/ORBdata.mat')
XJ2000 = mat["orb"]["XJ2000"]
T = mat["orb"]["t"]
XDF = pd.DataFrame(XJ2000,columns = ["X","Y","Z","VX","VY","VZ"])
sp.furnsh("input/LRO_ES_90_202003_GRGM900C_L600.bsp")
sp.furnsh("input/naif0012.tls")
start,stop = mat['orb']['epoch']['et']

### plot moon and traj
# fig1 = plt.figure()
# ax = plt.axes(projection='3d')

# #plot moon
# rM = 1738
# u = np.linspace(0, 2 * np.pi, 50)
# v = np.linspace(0, np.pi, 50)
# x = rM * np.outer(np.cos(u), np.sin(v))
# y = rM * np.outer(np.sin(u), np.sin(v))
# z = rM * np.outer(np.ones(np.size(u)), np.cos(v))
# ax.plot_surface(x, y, z,cmap = "gray",zorder=1)

# #plot traj
# l=None
# ax.plot(XJ2000[:l,0],XJ2000[:l,1],XJ2000[:l,2],"r",zorder=10,label="traj all")
# ax.set_xlabel('X (km)')
# ax.set_ylabel('Y (km)')
# ax.set_zlabel('Z (km)')

# # Plot initial position
# state = sp.spkezr("LRO",start,"J2000","NONE","MOON")[0]
# ax.plot(state[0],state[1],state[2],"b.",zorder=10,label="Initial pos")        #Starting point

# # Plot final position
# state = sp.spkezr("LRO",stop,"J2000","NONE","MOON")[0]
# ax.plot(state[0],state[1],state[2],"c.",zorder=10,label="Final Pos True")        #Ending point True

# # Plot final position Precise Prop
# EndP = XJ2000[-1]
# ax.plot(EndP[0],EndP[1],EndP[2],"y.",zorder=10,label="Final Pos All")        #Starting point All
# ax.legend()
# plt.axis("equal")
###

### Second graph
# acquire data
Tstep = mat['orb']['epoch']['span'][1] - mat['orb']['epoch']['span'][0]
Tspan = stop-start
Err = []
for _ in range(6):Err.append([])

for i in tqdm(range(len(XJ2000))):
    Time = T[i]
    state = sp.spkezr("LRO",Time,"J2000","NONE","MOON")[0]
    for j in range(6):
        Err[j].append(state[j]-XJ2000[i,j])

# plotting the 6 errors graphs
# fig2 = plt.figure(figsize=(18,4.8))
# axes = fig2.subplots(2,3)
# for i in range(6):
#     axes[i//3,i%3].plot(np.arange(0,Tspan+Tstep/2,Tstep),Err[i],
#                         label="RMSE = "+'{:.2e}'.format(np.sqrt((np.array(Err[i])**2).sum()/len(Err[i])))+[" km"," km/s"][i//3])
#     axes[i//3,i%3].set_xlabel('s')
#     axes[i//3,i%3].set_ylabel(f'error in {["X (km)","Y (km)","Z (km)","VX (km/s)","VY (km/s)","VZ (km/s)"][i]}')
#     axes[i//3,i%3].set_ylim([-.2,-.2e-3][i//3],[.2,.2e-3][i//3])
#     axes[i//3,i%3].plot([0,Tspan],[0,0])
    # axes[i//3,i%3].legend()
    
# plotting the 2 3D-errors graphs
fig2 = plt.figure()
axes = fig2.subplots(2,1)
Error3D = [np.sqrt((np.array(Err[:3])**2).sum(0))*1e3,np.sqrt((np.array(Err[3:])**2).sum(0))*1e5]
for i in range(2):
    axes[i].plot((T-T[0])/3600,Error3D[i],
                   label="RMSE = "+'{:.2f}'.format(np.sqrt((Error3D[i]**2).sum()/len(Error3D[i])))+[" m"," cm/s"][i])
    axes[i].set_xlabel('hours')
    axes[i].set_ylabel(f'error in {["position (m)","velocity (cm/s)"][i]}')
    axes[i].set_ylim(0,[300,30][i])
    axes[i].legend()

# close everything and plot
sp.kclear()
plt.show()

# sp.dafec(sp.spklef("input/LRO_ES_90_202003_GRGM900C_L600.bsp"),1000)