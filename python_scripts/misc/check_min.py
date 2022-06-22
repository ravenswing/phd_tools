import matplotlib.pyplot as plt
import numpy as np
import os
import sys


tol = 0.1 # percentage tolerance


def isMin (x,y,tol):
    ymax = max(y) 
    endi = len(data)-int(np.floor(len(data)/10)) 
    final10 = y[endi:]
    std = np.std(final10) 
    if std < np.absolute(ymax*tol*0.01):
        return True
    else:
        return False


if len(sys.argv) < 2:
    print("\nWhat file?\n")
    sys.exit()
with open(sys.argv[1]) as f:
    lines = f.readlines()
data = [l for l in lines if l[0] not in ("@", "#")]
data = [[float(val) for val in line.split()[:2]] for line in data]
x, y = [l[0] for l in data], [l[1] for l in data]
print( len(data) )
print( isMin(x,y,tol) )


'''
for line in lines:
    if line[0] == "@" and line.split()[1] == "title":
        plt.title(" ".join(line.split()[2:]).replace('"', ""))
    if line[0] == "@" and line.split()[1] == "xaxis":
        plt.xlabel(" ".join(line.split()[3:]).replace('"', ""))
    if line[0] == "@" and line.split()[1] == "yaxis":
        plt.ylabel(" ".join(line.split()[3:]).replace('"', ""))
plt.show()
'''
