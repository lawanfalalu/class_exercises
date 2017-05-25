#definition of a polynomial function
import numpy as n
def poly(min,max,step):
    coeff=[];f=[];guess=[];u=0
    degree = int(input("Enter the DEGREE : "))
    for i in range(degree+1):
        coeff.append(float(input("Enter coeffecient no. %d "%(i+1))))
    vrange = n.arange(min,max,step)
    for j in vrange:
        #for k in coeff
        sum=0
        for k in range(degree+1):
            sum = sum + coeff[k]*j**(degree-k)
        f.append(sum)
        if(u>1):
            if(f[u-1]*f[u])<=0:
                guess.append(j)
        u=u+1
        
        print("f(%f) = %f "%(j,sum))
    print(guess)

#Calling the method
poly(-10,10,1)