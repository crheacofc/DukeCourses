#Super basic making of structured mesh
import numpy as np

x = np.linspace(0,1,3)
y = np.linspace(0,1,3)

file = open('../outputs/basic.txt','w')
count = 0
for i in range(3):
    for j in range(3):
        if count != 0:
            file.write('\n')
        file.write(str(x[j]))
        file.write(" ")
        file.write(str(y[i]))
        count += 1
file.close()
