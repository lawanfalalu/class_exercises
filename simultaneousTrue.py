# -*- coding: utf-8 -*-
"""
Created on Mon May 29 12:46:43 2017

@author: FALALU
"""
#WORKING PERFECTLY CORRECT 
#simultaneous eqautions

def solve(matrix, mul=1):
    width = len(matrix)
    if width == 1:
        return mul * matrix[0][0]
    else:
        sign = -1
        total = 0
        for i in range(width):
            m = []
            for j in range(1, width):
                buff = []
                for k in range(width):
                    if k != i:
                        buff.append(matrix[j][k])
                m.append(buff)
            sign *= -1
            total += mul * solve(m, sign * matrix[0][i])
        return total



#end of my function
def simultaneous2():    
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
        r2 = float(input("constant value of equation %d : "%(i2+1)))
        re.append(r2)
    det = solve(coeff,1)
    print(coeff, " is a matrix and the determinant is :",det)
    if det == 0:
        print("no solution")
    else:
        for st in range(n):
            import copy
#            coeff2 = coeff[:]
#            coeff2=[]
            coeff2=copy.deepcopy(coeff) # nan wurin fa akwai wahala ashe-ashe baya copy
            print(coeff," is copied to",coeff2)
            for st2 in range(n):
                coeff2[st2][st]=re[st2]
            print("coeff undergoes change to :",coeff)
            print("X",st,"=",solve(coeff2,1)/det)
            
        




#    simultaneous(coeff,re)
simultaneous2()    
    
#print(coeff,re)