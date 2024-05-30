from pymatreader import read_mat
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spiceypy as sp
from datetime import datetime


mat = read_mat('C:/Users/quent/Documents/MATLAB/LHPOP/ORBdata.mat')
XJ2000 = mat["orb"]["XJ2000"]
XDF = pd.DataFrame(XJ2000,columns = ["X","Y","Z","VX","VY","VZ"])
sp.furnsh("lrorg_2021074_2021166_v01.bsp")
sp.furnsh("naif0012.tls")

### plot moon and traj
fig = plt.figure()
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
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
###

### plot altitude
# fig, ax = plt.subplots()
# plt.plot(np.sqrt(XDF["X"]**2+XDF["Y"]**2+XDF["Z"]**2)-1738)

### Plot initial position LRO 4/1/2021 12am
start = sp.str2et(str(datetime(2021,4,1)))
state = sp.spkezr("LRO",start,"J2000","NONE","MOON")[0]
ax.plot(state[0],state[1],state[2],"b.",zorder=10,label="Initial pos")        #Starting point

# Plot final position LRO 4/1/2021 1am
end = sp.str2et(str(datetime(2021,4,11)))
state = sp.spkezr("LRO",end,"J2000","NONE","MOON")[0]
ax.plot(state[0],state[1],state[2],"c.",zorder=10,label="Final Pos True")        #Ending point True

# Plot final position Fast Prop
EndF = np.array([ 7.63193760e+02,  1.16317361e+03, -1.15315866e+03, -1.18917603e+00,
       -3.35762839e-01, -1.10420559e+00])
ax.plot(EndF[0],EndF[1],EndF[2],"g.",zorder=10,label="Final Pos Fast")        #Starting point Fast

# Plot final position Precise Prop
EndP = XJ2000[-1]
ax.plot(EndP[0],EndP[1],EndP[2],"y.",zorder=10,label="Final Pos All")        #Starting point All

sp.kclear()
plt.legend()
plt.show()

def RMS(V1,V2):
    return np.sqrt(((V1-V2)[:3]**2).sum())

print("RMSF: ",RMS(state,EndF))
print("RMSP: ",RMS(state,EndP))