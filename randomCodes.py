# -*- coding: utf-8 -*-
"""
Created on Mon May 29 19:50:40 2017

@author: FALALU
"""

#x1 = (c- (A1*x2 +A2*x3))
def NewtonRaphson(coeff,dfdx,x,convergence,degree):
    k=0 
    if dfdx(x,coeff)==0:
        print("Halts!!! dfdx = 0")
        return
    print("X%d =%.5f , f(X%d)=%.5f , f'(X%d)=%.5f" % (k,x,k,f(x,degree,coeff),k,dfdx(x,coeff))) 
    while(abs(f(x,degree,coeff))>convergence): 
        k +=1 
        x=x-float(f(x,degree,coeff))/dfdx(x,coeff) 
        print("X%d =%.5f , f(X%d)=%.5f , f'(X%d)=%.5f" % (k,x,k,f(x,degree,coeff),k,dfdx(x,coeff))) 
         
    print ("The value of X=%.5f" %(x)) 
def f(X,degree,coeff): 
    Equation=0 
    indexofCoeff=0 
    for e in coeff: 
        Equation=Equation + e*pow(X,(degree-indexofCoeff)) 
        indexofCoeff +=1 
    return Equation 

def dfdx(X,coeff): 
    Equation2=0 
    indexofCoeff=0 
    for q in range(len(coeff),1,-1): 
        Equation2= Equation2 + (q-1)*coeff[indexofCoeff]*pow(X,q-2) 
        indexofCoeff +=1 
    return Equation2 
