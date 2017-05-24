#factoria fuction
def fact(n):
    fact = 1
    for i in range(1, n+1):
        fact = fact*i
    print("the factorial of %d is %.4E" % (n,fact))

#*************

#This is a function to solve a quadratic equation
import math
def QuadEq():
    print("Enter the coefficient of the Quadratic equation?")
    a = float(input("A = "))
    b = float(input("B = "))
    c = float(input("C = "))
    d =  pow(b,2)-(4*a*c)
    if(d<0):
        print("Complex root!!!")
    X1 = (-b + math.sqrt(d) )/(2*a)
    X2 = (-b + math.sqrt(d) )/(2*a)
    print("X1= %.4f and X2= %.4f" % (X1,X2) )
print("Good Bye!" )
QuadEq()

*****Simultaneous

def g(x):
    g = (-4-x**3 -3 * x**2)/6
    return(g)

x[0] = 0
for k in range(1000):
    x[k+1] = g(x[k])
    print(k+1, x[k+1])
    if(abs(x[x+1] - x[k]) < 0.1):
        break

*****another
def f(x):
    k = -x**3 + 3*x**2 + 6*x+4
    print(k)


for i in range(20):
    root =f(i)
    print(i,root)

U************
def f(x):
    k = -x**3 + 3*x**2 + 6*x+4
    print(x,k)
    
for i in range(20):
    root = f(i)
    if(i>0):
        product = f(i) * f(i-1)
        if(product<0):
            root = i-1
***************egg.py
from math import *
m = 47
M = 67
p = 1.038
c = 3.7
K = 5.4E-3
T0 = 4
Tw = 100
Ty = 70

t =  ( (pow(M,(2/3) * c * pow(p,(1/3)) )/( K * pow(pi,2) * pow( ((4*pi )/3), (2/3) ) ) ) * log( 0.76*( (T0-Tw)/(Ty-Tw) ) ) )
print("the time in seconds to boil the egg is: %.4E" % (t))