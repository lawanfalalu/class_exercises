# -*- coding: utf-8 -*-
"""
Created on Tue May 30 02:22:47 2017

@author: FALALU
"""

#Bisection method
#defining a function named Bisection
def Bisection(F,P1,P2,convergency,degree,coeff):
    aggregate = []
    Intermediate = []
    Xm_Next = []
#    cXmHolder = []
    j = 0
    withHeld1 = float(F(P1,degree,coeff))
    withHeld2 = float(F(P2,degree,coeff))
    
    #checking for sign changes
    if(withHeld1 * withHeld2)>0: 
       print("No solution eist for the inputed equation ")
    Xm = (P1 + P2) / 2.0
    Xu = F(Xm,degree,coeff)
    
    #checking condition for the convergence
    while(abs(Xu) > convergency):
        Xm_Next.append(Xu)
        if((withHeld1 * Xu)>0):
            withHeld1 = Xu
            P1 = Xm
            Intermediate.append(Xm)
        else:
            P2 = Xm
            aggregate.append(Xm)
        Xm=float(P1+P2)/2
        Xu=F(Xm,degree,coeff) 
    for e in  aggregate:
        print("The Values of F(X) within these points is: X%d=%.5f" % (j,e))
        j += 1
    print("The final root of the function is:  %.5f" %(e))
       
def F(X,degree,coeff):
    function = 0
    for e in range(degree + 1):
        function = function + coeff[e] * X**(degree-e)
    return function


def poly():
    guess = []
    f = []
    coeff = []
    degree = int(input("\nEnter the degree of the equation: "))
    for i in range(degree+1):
        k = 1 + i
        coeff.append(float(input("\nEnter the value of NO. %d coefficient of X: " % (k))))
    min = int(input("\nEnter Minimum range: "))
    max = int(input("\nEnter Maximum range: "))    
    
    u = 0
    for j in range(min,max):
        sum = 0
        for k in range(degree + 1):
            sum = sum + coeff[k] * j**(degree - k)
        f.append(sum)

        if(u>0):
            sign = (f[u-1]*f[u])
            if(sign <= 0):
                guess.append(j)
        u = u + 1
    
    print("\nThe guess points on X-axis is/are", guess)
    for x in range(len(guess)):
        X1 = int(input("\nChoose the two point values of X from this:  "+ str(guess)+":  "))
        X2 = int(input("\nChoose the two point values of X from this:  "+ str(guess)+":  ")) 
        convergency = float(input("\nEnter the maximum acceptable error: ")) 
    
        solution = Bisection(F,X1,X2,convergency,degree,coeff)
        print(solution)
poly()