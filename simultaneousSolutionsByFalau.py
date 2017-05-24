#simultaneous eqautions

import numpy as np
def simultaneous(p,q):
    A = np.array(p)
    B = np.array(q)
    X = np.linalg.solve(A,B)
    print("The solutions are:")
    for tmp in X:
        print(tmp)
#   print("X = ",X[0],"Y = ",X[1])
#end of my function

coeff=[]
n = int(input(" How many equations ? "))
for i in range(n):
    temp=[]
    print("getting coefficients of equation",i+1)
    for j in range(n):
        r = float(input("coefficient %d of equation %d : "%((j+1),(i+1))))
        temp.append(r)
    coeff.append(temp)
    temp=[]
re=[]
for i2 in range(n):
    r2 = float(input("resultant value of equation %d : "%(i2+1)))
    re.append(r2)
        
simultaneous(coeff,re)    
    
#print(coeff,re)