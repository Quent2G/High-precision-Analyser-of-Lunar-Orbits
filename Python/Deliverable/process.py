###Processing functions used in the other files

import numpy as np
from tqdm import tqdm
import spiceypy as sp
import pandas as pd

### Converts the spacecraft position from the inertial frame to the Earth-Moon rotational one
# Input:  S = Spacecraft states in inertial frame, output of propagator
#         T = Associated times, output of propagator
#         c = 1 if there is a converged trajectory (= optimization mode)
#         Sc= converged trajectory (in optimization mode) or not used
# Output: RS = Position of spacecraft in rotational frame
#         RE = Position of the Earth in rotational frame
#         RSC= Converged trajectory of spacecraft in rotational frame

def rotational(S,T,c,Sc):
    sp.furnsh("input/de430.bsp")
    RS,RE,RSC=[],[],[]
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

        RS.append(Mir@np.array(S[t,:3]))
        if c: RSC.append(Mir@np.array(Sc[t,:3]))
        RE.append([-norm(SE),0,0])
        
    RS = np.array(RS)
    RE = np.array(RE)
    if c: RSC = np.array(RSC)    
    else: RSC = 0
    sp.kclear()
    return RS,RE,RSC


### Get the position of the model CR3BP trajectory in the chosen frame
# Input:  path = Coordinates file of a CR3BP trajectory
#         RotationalF = 1 if the chosen frame is rotational (inertial otherwise)
#         T = Integration times, output of propagator
# Output: List of 3 lists: the list of X coordinates, same with Y and same with Z

def Model(path,RotationalF,T):
    # Read the CR3BP file and get the scaled coordinates
    LU = 389703
    TU = 382981
    mu = 1.215058560962404E-2
    with open(path) as file:
        Lines = [l.strip() for l in file.readlines()]

    Time,X,Y,Z = [],[],[],[]
    for i in range(1,len(Lines)):
        posRot = Lines[i].split(",")
        t,posRot = float(posRot[0])*TU,[float(x)*LU for x in posRot[1:4]]
        Time.append(t)
        X.append(posRot[0])
        Y.append(posRot[1])
        Z.append(posRot[2])
    
    # If the inertial frame was chosen, a fitting on ephemeris data has to be performed
    # followed by a conversion in the inertial and Moon-centered frame. 
    if not RotationalF:
        sp.furnsh("input/de430.bsp")
        RS = []
        i=0
        for t in Time:
            # Fitting first
            et = T[0]+t
            SE = sp.spkezr("EARTH",et,"J2000","NONE","MOON")[0]
            SB = (1-mu)*SE/norm(SE)*LU

            # Get the inertial frame vectors
            uR = -SE[:3]/norm(SE[:3])
            uTh = -SE[3:]/norm(SE[3:])
            uZ = np.cross(uR,uTh)

            # Conversion
            RS.append(SB[:3] + X[i]*uR + Y[i]*uTh + Z[i]*uZ)
            i+=1
        sp.kclear()
        return [r[0] for r in RS],[r[1] for r in RS],[r[2] for r in RS]

    # Ifelse the rotational frame was chosen, only a translation from barycentred frame to moon-centered frame is needed
    else: return X-(1-mu)*LU*np.ones(len(X)),Y,Z


### Allows to read a file taken from the Horizon website https://ssd.jpl.nasa.gov/horizons/app.html#/
def read_horizon(path):
    with open(path) as file:
        Lines = [l.strip() for l in file.readlines()]

    Data=[]

    #skip start
    l=1
    while Lines[l-1] != "$$SOE":
        l+=1

    while Lines[l] != "$$EOE":
        Data.append([])
        Data[-1].append(float(Lines[l].split()[0]))
        l+=1
        Pos = Lines[l].split("=")
        for i in range(3):
            Data[-1].append(float(Pos[i+1].split()[0]))
        l+=1
        Vel = Lines[l].split("=")
        for i in range(3):
            Data[-1].append(float(Vel[i+1].split()[0]))
        l+=1
        if "Orion." in path or "Orion2." in path:
            Vel = Lines[l].split("=")
            for i in range(3):
                Data[-1].append(float(Vel[i+1].split()[0]))
            l+=1
    

    col= ["T","X","Y","Z","VX","VY","VZ"]
    if "Orion." in path or "Orion2." in path: col += ["LT","RG","RR"]
    return(pd.DataFrame(Data,columns=col))

#Process the norm of a 6 item state
def norm(R):
    return np.sqrt((R**2).sum())

#Plot a state so that it is ready to use for Matlab LoadState.m function
def ConvMatlab(PosInit):
    print(f"""X = {PosInit[0]}; % X-Coordinate
        Y = {PosInit[1]};    % Y-Coordinate
        Z = {PosInit[2]};    % Z-Coordinate
        VX = {PosInit[3]};    % X-Velocity Coordinate
        VY = {PosInit[4]};    % Y-Velocity Coordinate
        VZ  = {PosInit[5]};    % Z-Velocity Coordinate""")
