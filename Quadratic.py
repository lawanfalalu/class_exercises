#Quadratic equautions
import math
def Quad():
    print("Getting the co-efficients...")
    A = float(input(" A = "))
    B = float(input(" B = "))
    C = float(input(" C = "))    
    D = B**2 -(4 * A * C)
    if (D < 0):
        print("Complex root!!!")
        return
    D = D**(0.5)
    X1 = (-B + D)/(2*A)
    X2 = (-B - D)/(2*A)
    print("X1 = %.4f and X2 = %.4f " %(X1,X2))
    print("GOOD BYE")