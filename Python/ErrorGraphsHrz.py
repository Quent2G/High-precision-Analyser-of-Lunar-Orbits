from pymatreader import read_mat
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spiceypy as sp
import datetime as dt
from tqdm import tqdm
from Read_Horizon import read_horizon

mat = read_mat('../MATLAB/LHPOP/output/ORBdata.mat')
XJ2000 = mat["orb"]["seq"]["a"]["XJ2000"]
T = mat["orb"]["seq"]["a"]["t"]
# XJ2000 = mat["orb"]["XJ2000"]
# T = mat["orb"]["t"]
XDF = pd.DataFrame(XJ2000,columns = ["X","Y","Z","VX","VY","VZ"])
start,stop = T[0],T[-1]
Hrz = read_horizon("input/Orion29_16.txt")
# Hrz = read_horizon("input/CAPSTONE_25_00__01_12.txt")
# Hrz = read_horizon("input/Change5.txt")

g=0
### plot moon and traj
if g:
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


    Hrz = read_horizon("input/OrionFull.txt")
    ax.plot(Hrz["X"],Hrz["Y"],Hrz["Z"],"r",zorder=10,label="Artemis I mission")
    Hrz = read_horizon("input/Orion29_16.txt")


    #plot traj
    # ax.plot(Hrz["X"],Hrz["Y"],Hrz["Z"],"r",zorder=10,label="Trajectory")
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')

    # Plot initial position
    state = Hrz.loc[0]
    ax.plot(state["X"],state["Y"],state["Z"],"b.",zorder=10,label="Initial pos Orion")        #Starting point

    # Plot final position
    state = Hrz.loc[len(Hrz)-1]
    ax.plot(state["X"],state["Y"],state["Z"],"c.",zorder=10,label="Final Pos Orion")        #Ending point True

    # Plot final position Precise Prop
    EndP = XJ2000[-1]
    ax.plot(EndP[0],EndP[1],EndP[2],"y.",zorder=10,label="Final Pos Propagation")        #Starting point All

    ax.legend()
    plt.axis('equal')
###

### Second graph
# acquire data
Tstep = T[1] - T[0]
Tspan = stop-start
Err = []
for _ in range(6):Err.append([])

for i in tqdm(range(len(XJ2000))):
    Time = T[i]
    state = np.array(Hrz.loc[i].drop("T"))
    for j in range(6):
        Err[j].append(state[j]-XJ2000[i,j])

# plotting the 2 3D-errors graphs
fig2 = plt.figure()
axes = fig2.subplots(2,1)
Error3D = [np.sqrt((np.array(Err[:3])**2).sum(0))*1e3,np.sqrt((np.array(Err[3:])**2).sum(0))*1e5]
for i in range(2):
    axes[i].plot((T-T[0])/86400,Error3D[i],
                   label="RMSE = "+'{:.2f}'.format(np.sqrt((Error3D[i]**2).sum()/len(Error3D[i])))+[" m"," cm/s"][i])
    axes[i].set_xlabel('days')
    axes[i].set_ylabel(f'error in {["position (m)","velocity (cm/s)"][i]}')
    # axes[i].set_ylim(0,[300,30][i])
    axes[i].legend()

# close everything and plot
plt.show()