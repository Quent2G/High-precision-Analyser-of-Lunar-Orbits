import pandas as pd

def read_horizon(path):
    with open(path) as file:
        Lines = [l.strip() for l in file.readlines()]

    Data=[]

    #skip debut
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

# read_horizon("input/CAPSTONE_25_00__01_12.txt")