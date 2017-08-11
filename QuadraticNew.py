# -*- coding: utf-8 -*-
"""
Created on Mon May 29 15:14:49 2017

@author: FALALU
"""

def Quadratic():
    import math
    print("------OF THE FORM Ax2 + Bx + C = 0 ------")
    a = float(input("Enter the value A :"))
    b= float(input("Enter the value of B :"))
    c = float(input("Enter the value of C :"))
    D = b**2 - 4*a*c
    if D < 0:
        rootD = math.sqrt(abs(D))
        r1 = complex(-b,rootD)/(2*a)
        r2 = complex(-b,-rootD)/(2*a)
    else:
        rootD = math.sqrt(D)
        r1 = (-b + rootD)/(2*a)
        r2 = (-b - rootD)/(2*a)
    print("The results are : X1 = ",r1," and X2 = ",r2)
Quadratic()