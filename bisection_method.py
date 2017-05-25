#Definition of Polynomial Function
index=-1;accum=[]
degree=0; coeff=[]; guess = [];
def Poly():
    u=0; sum2=0; F=[];global degree
    degree = int(input("Enter the degree:  "))
    for i in range (degree+1):
        c=float(input("Enter your coefficients  "))
        coeff.append(c)
    for j in range (-20, 20):
        sum1=0
        for k in range(degree+1):
            sum1=sum1+coeff[k]*j**(degree-k)
            
        print("f(%.0f)f=%.2f"%(j,sum1))
        F.append(sum1)
        if(u>0):
            sign=F[u-1]*F[u] 
            if(sign<=0):
                guess.append(j)
                
        u=u+1
    print("The Guess:" ,guess)

def func(je):
    fx=0
    #print(degree,coeff)
    for k in range(degree+1):
        fx=fx+coeff[k]*je**(degree-k)
    return fx
    
def bisect(a,b):
    global index
    global accum
    if index>=10:
        return
    avg=(a+b)/2
    accum.append(avg)
    index=index+1
    print("index is", index)
    print("cccc",accum)
    if index>0:
        if abs(accum[index]-accum[index-1])<0.0001:
            print("cccccc",accum[index])
    else:
        if (func(a)*func(avg)<=0):
            bisect(a,avg)
        else:
            bisect(avg,b)
    
        
Poly()
bisect(-1,0)