# coding: utf-8
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

data = np.loadtxt('data.txt')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(0,26):
	X = data[i][0]
	Y = data[i][1]
	Z = data[i][2]
	ax.scatter(X,Y,Z)



plt.show()
