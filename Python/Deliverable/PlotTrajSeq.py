import matplotlib.pyplot as plt
import numpy as np
from pymatreader import read_mat
import spiceypy as sp

fig = plt.figure(figsize=(10,8))
ax = plt.axes(projection='3d')
l=int(1e10)
# l=20000
c=0
model=0
earth=0

###Plot seq
mat = read_mat('../MATLAB/LHPOP/output/ORBdata.mat')

for e in list(mat["orb"]["seq"].keys())[1:]:
    XJ2000 = mat["orb"]["seq"][e]["XJ2000"]
    T = mat["orb"]["seq"][e]["t"]
    if c: XC = mat["orb"]["seq"][e]["XC"]
    ax.plot(XJ2000[:l,0],XJ2000[:l,1],XJ2000[:l,2],label="seq "+e)
    ax.scatter(XJ2000[0,0],XJ2000[0,1],XJ2000[0,2])
    if c: ax.plot(XC[:l,0],XC[:l,1],XC[:l,2],"orange",zorder=10,label="traj corrected")
    if c: ax.scatter(XC[0,0],XC[0,1],XC[0,2],c="orange")

ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')
ax.scatter(0,0,0,c="gray",label = "Moon")
plt.legend()
plt.axis("equal")
sp.kclear()