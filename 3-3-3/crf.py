import numpy as np
from itertools import product

x = np.zeros((3,3))
x[0][2] = 1
x[1][1] = 1
x[1][2] = 1
x[2][2] = 1
y = np.array([1,0,0])
z = np.ones((3,3))
z[0][0] = 0
z[2][0] = 0
z[2][1] = 0

theta = 0.8
mu = 1.2

num = theta**7 * mu

denum = 0.0
for matrix in product([0, 1], repeat=9):
	config = np.array([matrix[0:3],matrix[3:6],matrix[6:9]])
	tmp = 1
	for i in range(3):
		for j in range(3):
			if config[i][j] == z[i][j]:
				tmp = tmp * theta
			if config[i][j] == 1 and y[i] == 1:
				tmp = tmp * mu
	denum += tmp

print num, denum, num/denum


	