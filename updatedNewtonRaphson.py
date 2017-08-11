# -*- coding: utf-8 -*-
"""
Created on Mon May 29 17:44:44 2017

@author: FALALU
"""
degree =1
def Raphson():
    degree = 2
    coeff =[1,0,-3]
    
    
def NewRaph(f,dfdx,x,convergence,degree): 
    k=0 
    print("Intermediate Values of X are: X%d =%.5f , f(X%d)=%.5f , f'(X%d)=%.5f" % (k,x,k,f(x,degree),k,dfdx(x))) 
    while(abs(f(x,degree))>convergence): 
        k +=1 
        x=x-float(f(x,degree))/dfdx(x) 
        print("Intermediate Values of X are: X%d =%.5f , f(X%d)=%.5f , f'(X%d)=%.5f" % (k,x,k,f(x,degree),k,dfdx(x))) 
         
    print ("The value of X=%.5f" %(x)) 
def f(X,degree): 
    Equation=0 
    indexofCoeff=0 
    for e in coeff: 
        Equation=Equation + e*pow(X,(degree-indexofCoeff)) 
        indexofCoeff +=1 
    return Equation 

def dfdx(X): 
    Equation2=0 
    indexofCoeff=0 
    for q in range(len(coeff),1,-1): 
        Equation2= Equation2 + (q-1)*coeff[indexofCoeff]*pow(X,q-2) 
        indexofCoeff +=1 
    return Equation2 
 
         
NewRaph(f,dfdx,20,0.001,degree) 