from math import pow,pi,log
M = 67
p = 1.038
c = 3.7
K = 5.4E-3
T0 = 4
Tw = 100
Ty = 70
num = pow(M,(2/3))*c*pow(p,(1/3))
den = K*pow(pi,2)*pow(((4*pi)/3),(2/3))
logVal = log( 0.76*((T0-Tw)/(Ty-Tw))) 
t = (num/den)*logVal
print("It takes",round(t,2),"to boil the EGG")