import matplotlib.pyplot as plt
import numpy as np
from pymatreader import read_mat
import spiceypy as sp

fig = plt.figure(figsize=(10,8))
ax = plt.axes(projection='3d')
l=int(1e10)
# l=20000

###Plot DRO_NRHOSimuated
mat = read_mat('../MATLAB/LHPOP/output/ORBdata.mat')
XJ2000 = mat["orb"]["XJ2000"]
T = mat["orb"]["t"]
ax.plot(XJ2000[:l,0],XJ2000[:l,1],XJ2000[:l,2],"r",zorder=10,label="traj sim")
ax.scatter(XJ2000[0,0],XJ2000[0,1],XJ2000[0,2],c="r")
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')
ax.scatter(0,0,0,c="gray")

### Plot DRO_NRHOTheoric
sp.furnsh("input/de430.bsp")
et26_11 = 722692869.182957
SE = sp.spkezr("EARTH",et26_11,"J2000","NONE","MOON")[0]
def norm(R):
    return np.sqrt((R**2).sum())
LU = 389703
TU = 382981
mu = 1.215058560962404E-2
path = "input/statesNRHOCapstone.csv"
with open(path) as file:
    Lines = [l.strip() for l in file.readlines()]

T,X,Y,Z = [],[],[],[]
for i in range(1,len(Lines)):
    posRot = Lines[i].split(",")
    t,posRot = float(posRot[0])*TU,[float(x)*LU for x in posRot[1:4]]
    T.append(t)
    X.append(posRot[0])
    Y.append(posRot[1])
    Z.append(posRot[2])

def norm(R):
    return np.sqrt((R**2).sum())

# Convert position from rotational frame to inertial at each step
i=0
RS = []
RE = []
REM = []
for j in range(1):
    for t in T[:int(l*len(T)/len(XJ2000))]:
        et = 722692869.182957+t+j*T[-1]
        SE = sp.spkezr("EARTH",et,"J2000","NONE","MOON")[0]
        SB = (1-mu)*SE/norm(SE)*LU

        uR = -SE[:3]/norm(SE[:3])
        uTh = -SE[3:]/norm(SE[3:])
        uZ = np.cross(uR,uTh)

        RE.append(SE[:3])
        REM.append(-LU*uR)
        RS.append(SB[:3] + X[i]*uR + Y[i]*uTh + Z[i]*uZ)
        i+=1
    i=0           #For second period
# for t in T[:int(l*len(T)/len(XJ2000))]:
#     et = 722692869.182957+t+T[-1]
#     SE = sp.spkezr("EARTH",et,"J2000","NONE","MOON")[0]
#     SB = (1-mu)*SE/norm(SE)*LU

#     uR = -SE[:3]/norm(SE[:3])
#     uTh = -SE[3:]/norm(SE[3:])
#     uZ = np.cross(uR,uTh)

#     RE.append(SE[:3])
#     REM.append(SB[:3]-mu*LU*uR)
#     RS.append(SB[:3] + X[i]*uR + Y[i]*uTh + Z[i]*uZ)
#     i+=1

RS = np.array(RS)
RE = np.array(RE)
REM = np.array(REM)
ax.plot(RS[:,0],RS[:,1],RS[:,2],c="g",label="satellite")
ax.scatter(RS[0,0],RS[0,1],RS[0,2],c="g")
ax.plot(RE[:,0],RE[:,1],RE[:,2],c="b",label="Earth")
ax.scatter(RE[0,0],RE[0,1],RE[0,2],c="b")
ax.plot(REM[:,0],REM[:,1],REM[:,2],c="c",label="Earth_CR3BP")


# theta = np.linspace(0, 2 * np.pi, 100)
# x = np.cos(theta)*LU
# y = np.sin(theta)*LU
# ax.plot(x,y)

plt.legend()
plt.axis("equal")
sp.kclear()