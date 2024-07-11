import matplotlib.pyplot as plt
import numpy as np
from pymatreader import read_mat
import spiceypy as sp
from tqdm import tqdm

fig = plt.figure(figsize=(10,8))
ax = plt.axes(projection='3d')
l=int(1e10)
# l=20000
c=0
model=0
earth=0

###Plot DRO_NRHOSimuated
mat = read_mat('../MATLAB/LHPOP/output/ORBdata.mat')
XJ2000 = mat["orb"]["XJ2000"]
T = mat["orb"]["t"]
sp.furnsh("input/de430.bsp")
if c: XC = mat["orb"]["XC"]

# Convert position from inertial frame to rotational at each step
RS = []
if c: RSC = []
RE = []
LU = 389703
TU = 382981
mu = 1.215058560962404E-2
et26_11 = 722692869.182957

def norm(R):
    return np.sqrt((R**2).sum())

for t in tqdm(range(len(T))):
    et = T[t]
    SE = sp.spkezr("EARTH",et,"J2000","NONE","MOON")[0]

    uR = -SE[:3]/norm(SE[:3])
    uTh = -SE[3:]/norm(SE[3:])
    uZ = np.cross(uR,uTh)

    Mri = np.zeros((3,3))
    Mri[:,0] = uR
    Mri[:,1] = uTh
    Mri[:,2] = uZ
    Mir = np.linalg.inv(Mri)

    RS.append(Mir@np.array(XJ2000[t,:3]))
    if c: RSC.append(Mir@np.array(XC[t,:3]))
    RE.append([-norm(SE),0,0])

RS = np.array(RS)
if c: RSC = np.array(RSC)
RE = np.array(RE)
# ax.plot(RS[:,0],RS[:,1],RS[:,2],"g",zorder=10,label="NRHO - Gateway")
ax.plot(RS[:,0],RS[:,1],RS[:,2],"r",zorder=10,label="DRO - Orion")
ax.scatter(RS[0,0],RS[0,1],RS[0,2],c="r")
if c: ax.plot(RSC[:,0],RSC[:,1],RSC[:,2],"orange",zorder=10,label="traj corrected")
if c: ax.scatter(RSC[0,0],RSC[0,1],RSC[0,2],c="orange")
if earth: ax.scatter(-LU,0,0,c="b",label="Earth_CR3BP")
if earth: ax.plot(RE[:,0],RE[:,1],RE[:,2],c="c",label="Earth")
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')
ax.scatter(0,0,0,c="gray",label = "Moon")
plt.axis("equal")

### Plot DRO_NRHOTheoric
LU = 389703
TU = 382981
mu = 1.215058560962404E-2
# path = "input/statesDRO14.csv"
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

if model: ax.plot(X-(1-mu)*LU*np.ones(len(X)),Y,Z,c="g",label="CR3BP")
if model: ax.scatter(X[0]-(1-mu)*LU,Y[0],Z[0],c="g")

plt.legend()
plt.axis("equal")
sp.kclear()