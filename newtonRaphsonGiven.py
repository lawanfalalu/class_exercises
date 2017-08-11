# -*- coding: utf-8 -*-
"""
Created on Mon May 29 03:36:34 2017

@author: FALALU
"""

# -*- coding: utf-8 -*- 
""" 
Created on Wed May 24 15:01:05 2017 
 
@author: user 
""" 
 
#newton Raphson method 
#coeff=[] 
#degree=int(input("Enter the degree of the quadratic equation to be solved: ")) 
#for i in range(degree+1): 
#    k=1+i 
#    c=float(input("Enter the value of %dth coefficient of X: " % (k))) 
#    coeff.append(c) 
#    
#def NewtonRaphson():
#    print("Calling newtonRaphson")
    
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
 
def Bisect(Xl,Xu,coeff,Xm):
    degree = len(coeff)-1
#    print(middle)
    middle.append(Xm)
    data = input("press any key to continue...")
    if len(middle)>1:
        for k in range(len(middle)):
            if abs(middle[k-1])-abs(middle[k]) < 0.0000001:
                print(middle[k],"is the solution")
                return
    print(f(Xl,degree,coeff)*f(Xu,degree,coeff))
    if (f(Xl,degree,coeff)*f(Xu,degree,coeff)) < 0:
        Xm = (Xl+Xu)/2
        middle.append(Xm)
        if (f(Xl,degree,coeff)*f(Xm,degree,coeff))<0:
            Xl = Xl
            Xu = Xm
        elif (f(Xl,degree,coeff)*f(Xm,degree,coeff))>0:
            Xl = Xm
            Xu = Xu
  
    Bisect(Xl,Xu,coeff,Xm)
#            print("the root is : ",Xm)
    
    
Bisect(-1,0,[1,-3,-4],0)
            
    
        
        
    
NewtonRaphson([1,0,-3],dfdx,20,1e-5,degree) 