#Gauss Iteration(Gauss-Siedel)
print("j\tx1\tx2\tx3")
x1=4; x2=-8;x3=15
x11=[]; x22=[];x33=[]
for j in range(500):
    x1=(4+x2-3*x3)/2
    x11.append(x1)
    x2=(-8-x1+2*x3)/9
    x22.append(x2)
    x3=(15-4*x1+8*x2)/10
    x33.append(x3)
    print(j,"\t", x1,"\t", x2, "\t", x3)
    if(j>0):
        if(abs(x11[j]-x11[j-1])<0.0000000000000001) and (abs(x22[j]-x22[j-1])<0.0000000000000001) and (abs(x33[j]-x33[j-1])<0.0000000000000001):
            print("x1, x2, x3 converged at:",x1,x2,x3)
            break