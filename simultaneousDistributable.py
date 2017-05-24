#simultaneous equations

import numpy as n

A = n.array([ [3,-4], [2,-3] ])
B = n.array([16,11])
X = n.linalg.solve(A,B)
print("X = ",X[0],"Y = ",X[1])