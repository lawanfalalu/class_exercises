from math import *
M = 67
p = 1.038
c = 3.7
K = 5.4E-3
T0 = 4
Tw = 100
Ty = 70
t = ((pow(M,(2/3))*c*pow(p,(1/3)))/(K*pow(pi,2)*pow(((4*pi)/3),(2/3))))  *  log( 0.76*((T0-Tw)/(Ty-Tw))) 

print("the time in seconds to boil the egg is: %.3E" % (t))