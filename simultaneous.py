#simultaneous eqautions

import numpy as np
def simultaneous(cXeq1,cYeq1,valEq1,cXeq2,cYeq2,valEq2):
    A = np.array([ [cXeq1,cYeq1], [cXeq2,cYeq2] ])
    B = np.array([valEq1,valEq2])
    X = np.linalg.solve(A,B)
    print("X = ",X[0],"Y = ",X[1])
    
# Getting the values
def getters():
    v1 = float(input("Enter the 1st coefficient in Eqn1 :"))
    v2 = float(input("Enter the 2nd coefficient in Eqn1 :"))
    v3 = float(input("Enter the value in Eqn1 :"))
    v4 = float(input("Enter the 1st coefficient in Eqn2 :"))
    v5 = float(input("Enter the 2nd coefficient in Eqn2 :"))
    v6 = float(input("Enter the value in Eqn2 :"))
    simultaneous(v1,v2,v3,v4,v5,v6)
    
getters()