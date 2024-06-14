import matplotlib.pyplot as plt

with open("input/Study.txt") as file:
    Lines = file.readlines()

fig,axS = plt.subplots()
axS.set_title(Lines[0][7:-1])
axS.set_xlabel(Lines[1][8:-1])
axS.set_ylabel(Lines[2][8:-1])

n = len(Lines)-5
X,Y = [],[]
for i in range(n):
    data = [float(x) for x in Lines[i+5].split()]
    X.append(data[0])
    Y.append(data[1])
axS.plot(X,Y)