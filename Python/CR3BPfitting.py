# Converts the initial state from a CR3BP orbit chosen on https://ssd.jpl.nasa.gov/tools/periodic_orbits.html#/intro
# in the Earth-Moon rotationnal and barycentered frame into the ephemeris Moon centered inertial frame 
# in order to use it as an initial state for the propagator.
#
# Input:  path = A CR3BP orbit csv file
# Output: The converted initial state ready to paste into the propagator LoadState.m function
#####

import matplotlib.pyplot as plt
import spiceypy as sp
import numpy as np
from process import norm, ConvMatlab

### Get original state
# path = "input/statesDRO14.csv"
path = "input/statesNRHOCapstone.csv"

### Read initial state
with open(path) as file:
    line = file.readline()
    line = file.readline()
X,Y,Z,VX,VY,VZ = [float(x) for x in line.strip().split(",")[1:]]

### Get Moon-Barycenter (MB) vector (CR3BP equivalent)
sp.furnsh("input/de430.bsp")
sp.furnsh("input/naif0012.tls")
et = sp.datetime2et(sp.datetime(2024,1,1))
ME = sp.spkezr("EARTH",et,"J2000","NONE","MOON")[0]
LU = 389703
TU = 382981
mu = 1.215058560962404E-2
MB = (1-mu)*ME/norm(ME)*LU

### Get rotational frame and retrieve the file coordinates
uR = -ME[:3]/norm(ME[:3])
uTh = -ME[3:]/norm(ME[3:])
uZ = np.cross(uR,uTh)
X = np.array(X)*LU
Y = np.array(Y)*LU
Z = np.array(Z)*LU
VX = np.array(VX)*LU/TU
VY = np.array(VY)*LU/TU
VZ = np.array(VZ)*LU/TU

### Calc Omg=dur/(dt*uTh)
dt = 1
MEm = sp.spkezr("EARTH",et-dt,"J2000","NONE","MOON")[0]
uRm = -MEm[:3]/norm(MEm[:3])
Omg = norm(uR-uRm)/dt/norm(uTh)

### Convert state from rotational barycentered frame to Moon-centered inertial
BS = X*uR + Y*uTh + Z*uZ
RS = MB[:3] + BS
VS = MB[3:] + VX*uR + VY*uTh + VZ*uZ + np.cross(Omg*uZ,BS)
init = list(RS)+list(VS)
ConvMatlab(init)
sp.kclear()