import matplotlib.pyplot as plt
import spiceypy as sp
import numpy as np
from ConvMatlab import ConvMatlab

# Get original state
path = "input/statesNRHOCapstone.csv"
with open(path) as file:
    line = file.readline()
    line = file.readline()
X,Y,Z,VX,VY,VZ = [float(x) for x in line.strip().split(",")[1:]]
# state0 = [float(x) for x in line.strip().split(",")[1:]]

def norm(R):
    return np.sqrt((R**2).sum())

sp.furnsh("input/de430.bsp")
et26_11 = 722692869.182957
ME = sp.spkezr("EARTH",et26_11,"J2000","NONE","MOON")[0]
LU = 389703
TU = 382981
mu = 1.215058560962404E-2
MB = (1-mu)*ME/norm(ME)*LU

uR = -ME[:3]/norm(ME[:3])
uTh = -ME[3:]/norm(ME[3:])
uZ = np.cross(uR,uTh)
X = np.array(X)*LU
Y = np.array(Y)*LU
Z = np.array(Z)*LU
VX = np.array(VX)*LU/TU
VY = np.array(VY)*LU/TU
VZ = np.array(VZ)*LU/TU
#Modify the VYo = -0.1031873007851023
# VY = -0.1
# VX = -0.01
# Z += 500

# Calc Omg=dur/(dt*uTh)
dt = 1
MEm = sp.spkezr("EARTH",et26_11-dt,"J2000","NONE","MOON")[0]
uRm = -MEm[:3]/norm(MEm[:3])
Omg = norm(uR-uRm)/dt/norm(uTh)

# Convert position from rotational frame to inertial
BS = X*uR + Y*uTh + Z*uZ
RS = MB[:3] + BS
VS = MB[3:] + VX*uR + VY*uTh + VZ*uZ + np.cross(Omg*uZ,BS)
init = list(RS)+list(VS)
ConvMatlab(init)
sp.kclear()
