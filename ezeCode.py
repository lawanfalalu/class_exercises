#definition of a polynomial function
import numpy as n
coeff=[];f=[];guess=[]; degree=0
f11=0
f22=0
a=0
u=0
ans = []
def poly(min,max,step):
    u=0
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
            #if((f[u-2]*f[u])<0 and (f[u-1]*f[u])!=0) :
            if(f[u-1]*f[u])<0 or f[u]==0:
                guess.append(j-1)
        u=u+1
        
        print("f(%f) = %f "%(j,sum))
    print(guess)

x1=1
newX1=0

def getf11(x1):
    x1 = guess[0];
    ans.append(x1)
    f11=0
    for k in range(degree+1):
        f11= f11 + coeff[k]*x1**(degree-k)
    return f11

def getf22(x1):
    f11=1
    for k in range(degree):
        f11= f11 + (degree-k)*coeff[k]*x1**(degree-k-1)
    return f11

       
def solvPoly():
    
    x1=guess[u]
    a=0 
    newX1 = x1-getf11(x1) / getf22(x1)
    ans.append(newX1)
    x1=newX1
    newX1 = x1-getf11(x1) / getf22(x1)
    ans.append(newX1)
    
    while(a>0 and (abs(ans[a]-ans[a-1]))>0.0001):
        x1=newX1
        a=a+1
        newX1 = x1-getf11(x1) / getf22(x1)
        ans.append(newX1)
        print("converged at",newX1)
    a=a+1
    
                
        
            
        #print(ans)
            
        
    print(ans)
#Calling the method
poly(-10,10,0.1)
solvPoly()