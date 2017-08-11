#Gauss Iteration(Gauss-Siedel)
print("J\tX1\tX2\tX3")
x1=4; x2=-8;x3=15
Y1=[]; Y2=[];Y3=[]
for j in range(500):
    x1=(4+x2-3*x3)/2
    Y1.append(x1)
    x2=(-8-x1+2*x3)/9
    Y2.append(x2)
    x3=(15-4*x1+8*x2)/10
    Y3.append(x3)
    print(j,"\t", x1,"\t", x2, "\t", x3)
    if(j>0):
        if(abs(Y1[j]-Y1[j-1])<1e-16) and (abs(Y2[j]-Y2[j-1])<1e-16) and (abs(Y3[j]-Y3[j-1])<1e-16):
            print("x1, x2, x3 converged at:",x1,x2,x3)
            break