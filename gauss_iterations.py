#Gauss iterartions
# In Gauss Seidel you don't need to define X1
# But in Gauss Jacobi X1 was defined
print("j \t X1 \t X2 \t X3")
y1=[];y2=[];y3=[]
x1=4
x2=-8
x3=15
for j in range(50):
    x1 = (4+x2-3*x3)/2
    y1.append(x1)
    x2 = (-8-x1+2*x3)/9
    y2.append(x2)
    X3 = (15 - 4*x1 + 8*x2)/10
    y3.append(x3)
    print(j,"\t",x1,"\t",x2,"\t",x3)
    if(j>0):
        if((abs(y1[j]-y1[j-1]))<0.00001) and ((abs(y2[j]-y2[j-1]))<0.00001) and ((abs(y3[j]-y3[j-1]))<0.00001):
            print("Converges at ",j)
            break