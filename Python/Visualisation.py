# Allows to plot a propagation after processing, the user can choose the frame, 
# or choose to add the converged trajectory if there is one, same with a CR3BP trajectory.
#
# Input:  mat = A propagation file (matlab file, just need to run HALO)
#       & SConv = A converged propagation file (optional, for optimization mode)
#       & path = A CR3BP trajectory (optional, for optimization mode)
# Output:
#         Plot the trajectory (user's choice)
#       & Plot the error graphs between the propagation and the real evolution
#####

import matplotlib.pyplot as plt
import numpy as np
from pymatreader import read_mat
import spiceypy as sp
from process import rotational,Model

fig = plt.figure(figsize=(10,8))
ax = plt.axes(projection='3d')

### User's choice
RotationalF = 0                                 #Plot the sequential in the Earth-Moon rotational frame (0=inertial)
Converged =   0                                 #Plot the initial and converged trajectory (for optimization mode)
earth =       0                                 #Plot the earth (only in rotat. frame)
model =       0                                 #Plot the CR3BP trajectory (for optimization mode)
# path = "input/statesDRO14.csv"                  #Path of the model (initial condition for optimization)
path = "input/statesNRHOCapstone.csv"

### Load the computed data
mat = read_mat('../MATLAB/HALO/output/ORBdata.mat')

for e in list(mat["orb"]["seq"].keys())[1:]:
    Ssat = mat["orb"]["seq"][e]["XJ2000"]               #States of the satellite in J2000
    T = mat["orb"]["seq"][e]["t"]                       #Time during sequence
    if Converged: SConv = mat["orb"]["seq"][e]["XC"]    #States of the converged trajectory
    else: SConv = 0

    ###Process the change of frame and the model if needed
    if RotationalF: Ssat,SEarth,SConv = rotational(Ssat,T,Converged,SConv)
    if model: ModP = Model(path,RotationalF,T)
    Plot = Ssat

    ###Plot the propagation
    ax.plot(Plot[:,0],Plot[:,1],Plot[:,2],label="seq "+e)
    ax.scatter(Plot[0,0],Plot[0,1],Plot[0,2])

    ###Plot depending on user's choice
    if Converged: 
        ax.plot(SConv[:,0],SConv[:,1],SConv[:,2],"orange",zorder=10,label="Traj. converged")
        ax.scatter(SConv[0,0],SConv[0,1],SConv[0,2],c="orange")
    if model: 
        ax.plot(ModP[0],ModP[1],ModP[2],c="g",label="CR3BP")
        ax.scatter(ModP[0][0],ModP[1][0],ModP[2][0],c="g")
    if earth:
        if not RotationalF:
            sp.furnsh("input/de430.bsp")
            SEarth = np.zeros((len(T),3))
            for i in range(len(T)):
                SEarth[i] = sp.spkezr("EARTH",T[i],"J2000","NONE","MOON")[0][:3]
        else: ax.scatter(-389703,0,0,c="b",label="Earth_CR3BP")
        ax.plot(SEarth[:,0],SEarth[:,1],SEarth[:,2],c="c")
        ax.scatter(SEarth[0,0],SEarth[0,1],SEarth[0,2],c="c")
if earth: ax.plot(SEarth[:,0],SEarth[:,1],SEarth[:,2],c="c",label="Earth")
        
###Graph
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')
ax.scatter(0,0,0,c="gray",label = "Moon")
plt.legend()
plt.axis("equal")
sp.kclear()
plt.show()