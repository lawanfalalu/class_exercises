# -*- coding: utf-8 -*-
"""
THIS IS ASSIGNMENT ONE OF WHICH IT CONTAINS SOLUTIONS TO:
    PERFECT EGG
    LINEAR EQUATIONS
        SIMULTANEOUS
            CRAMMERS RULE
            GAUSS SEIDEL ITERATIONS
            GAUSS JACOBI ITERATIONS
    NON-LINEAR EQUATIONS
        QUADRATIC
        POLYNOMIAL
            NEWTON RAPSON
            BISECTION METHOD
    AUST GRADING POINT
    LU DECOMPOSITIONS
    
Created on Sat May 27 17:12:48 2017

@author: FALALU
"""
###################################### AUST GRADING SYSTEM #############################################
import math
def GradePoints(mark):
    gSys = {"A":4.0,"A-":3.75,"B+":3.25,"B":3.0,"B-":2.75,"C+":2.25,"C":2.0,"C-":1.75,"D":1.0,"F":0.0,}
#    print(gSys[mark])
    return gSys[mark]
def GradeMe(mks):
    if mks>=95 and mks<=100:
        return "A"
    elif mks>=89:
        return "A-"
    elif mks>=83:
        return "B+"
    elif mks>=77:
        return "B"
    elif mks>=71:
        return "B-"
    elif mks>=65:
        return "C+"
    elif mks>=59:
        return "C"
    elif mks>=53:
        return "C-"
    elif mks>=48:
        return "D"
    else:
        return "F"
   
def GetGrade():
    trans={}
    points=0
#    credit= float(input("WHAT'S THE CREDIT UNIT ? :"))
#    print("\n")
    total= int(input("HOW MANY COURSES ? :"))
    if total > 0:
        for c in range(total):
            code = input("COURSES CODE %d ? :"%(c+1))
            mark = float(input("MARKS OPTAINED IN \"%s\" ? :"%(code)))
            trans[code]=mark

        for key,value in trans.items():
            print(key,value,GradeMe(value))
            points = points + GradePoints(GradeMe(value))
            print("\tCGPA = #",round((points/total),3))
    else:
        print("CANNOT CALCULATE %d COURSES !"%(total))
######################################END OF AUST GRADING SYSTEM #############################################

######################################BOIL PERFECT EGG #############################################
def PerfectEgg():  
    from math import pow,pi,log
    M = 67
    p = 1.038
    c = 3.7
    K = 5.4E-3
#    T0 = 4
    Tw = 100
    Ty = 70
    T0 = float(input("What's the Original Tempreture, T0 ?:"))
    num = pow(M,(2/3))*c*pow(p,(1/3))
    den = K*pow(pi,2)*pow(((4*pi)/3),(2/3))
    logVal = log( 0.76*((T0-Tw)/(Ty-Tw))) 
    t = (num/den)*logVal
    print("It takes",round(t,2)," MIN(s) to boil the EGG\n")

######################################END OF BOIL PERFECT EGG #############################################

###################################### QUADRATIC EQUATIONS #############################################


def Quadratic():             
    print("------OF THE FORM Ax2 + Bx + C = 0 ------")
    a = float(input("Enter the value A :"))
    b= float(input("Enter the value of B :"))
    c = float(input("Enter the value of C :"))
    D = b**2 - 4*a*c
    if D < 0:
        rootD = (abs(D))**0.5
        r1 = complex(-b,rootD)/(2*a)
        r2 = complex(-b,-rootD)/(2*a)
        print("The results are : X1 = ",r1," and X2 = ",r2)
    elif D == 0:
        print("The result is : X1 & X2 both = ",(-b/(2*a)))
    else:
        rootD = D**0.5
        r1 = (-b + rootD)/(2*a)
        r2 = (-b - rootD)/(2*a)
        print("The results are : X1 = ",r1," and X2 = ",r2)
    
    print("Press 1 to solve new Qaudratic Equation")
    print("Press 0 to GO TO the previous MENU")
    CHOICE=int(input("PLEASE ENTER YOUR CHOICE [0...1]:"))
    if CHOICE == 1:
        Quadratic()

    while CHOICE != 0:
         print("Press 1 to solve new Qaudratic Equation")
         print("Press 0 to GO TO the previous MENU")
         CHOICE=int(input("PLEASE ENTER YOUR CHOICE [0...1]:"))
         if CHOICE == 0:
             break
         elif CHOICE == 1:
             Quadratic()
        

    

        


###################################### END QUADRATIC EQUATIONS #############################################

####################################### Gauss Iteration(Gauss-Siedel) ######################################

                    
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 23:21:57 2017

@author: FALALU
"""

from math import *
m = []
def ML(m):
    if(isSquare(m)):
        zeros_diagonal = zeros_matrix(len(m), len(m))        
        for i in range(1, len(m), 1):
            for j in range(0, i):
                zeros_diagonal[i][j] = - m[i][j]
        return zeros_diagonal
    else:
        print ("Must be a square matrix")

def MU(m):
    if(isSquare(m)):
        zeros_diagonal = zeros_matrix(len(m), len(m))        
        for i in range(0, len(m)):
            for j in range(i+1, len(m)):
                zeros_diagonal[i][j]= -m[i][j]
        return zeros_diagonal
    else:
        print ("Must be a square matrix")
        
def diagonial_matrix(m):
    if(isSquare(m)):
        zeros_diagonal = zeros_matrix(len(m), len(m))
        for i in range(len(m)):
            for j in range(len(m)):
                zeros_diagonal[i][i] = m[i][i]
        return zeros_diagonal
    else:
        print ("Must be a square matrix")


def vector_vector_mult(v1, v2):
    result = 0
    for i in range(len(v1)):
        result = result + v1[i]*v2[i]
    return result

def V2Vadd(v1, v2):
    result = zeros_array(len(v1))
    for i in range(len(v1)):
        result[i] = v1[i]+v2[i]
    return result

def M_inverse(m):
    return matrix_scalar_mult(matrix_adjunct(m), (1/matrix_determinant(m)))

def matrix_scalar_mult(m, s):    
    for i in range(len(m)):
        for j in range(len(m[0])):
            m[i][j] = m[i][j] * s
    return m

def matrix_cofactor(m):
    result = []
    for i in range(len(m)):
        aRow=[]
        for j in range(len(m)):
            x=((-1)**(i+j))*matrix_determinant(adjunct_sub_matrix(m,i,j))
            aRow.append(x)
        result.append(aRow)
    return result

def isSquare(m):
    return  all(len(row) == len(m) for row in m)

def isList(a):
    return isinstance(a,list)

def dot_product(a1, a2):
    if (isList(a1) and isList(a2)):
        if (len(a1)==len(a2)):
            sum = 0
            for i in range(len(a1)):
                sum = sum + a1[i]*a2[i]
            return sum
        else:
            raise ValueError('a1 and a2 must be arrays of  same length')            
    else:
        raise ValueError('a1 and a2 must be must lists')   
    
def MatrixMult(m, v):
    result = []
    for row in m:
        aProduct = dot_product(row, v)
        result.append(aProduct)
    return result

def matrix_transpose(m):
    return [[x[i] for x in m] for i in range(len(m[0]))]

def matrix_adjunct(m):
    if(isSquare(m)):
        return matrix_transpose(matrix_cofactor(m))
    else:
        raise ValueError('Must be a square matrix')           

def matrix_determinant(m):
    return sum(recursion(m))


def zeros_array(size):
    result=[]
    for i in range(size):
        result.append(0)
    return result           

def zeros_matrix(dim1, dim2):
    result=[]
    for i in range(dim1):
        aRow = []
        for j in range(dim2):
            aRow.append(0)
        result.append(aRow)
    return result  

def adjunct_sub_matrix(m, exclu_row, exclu_col):
    sub_matrix=[]
    for i in range(len(m)):
        if(i == exclu_row):
            continue
        else:
            aRow =[]
            for j in range(len(m[i])):
                if(j == exclu_col):
                    continue
                else:
                    aRow.append(m[i][j])
            sub_matrix.append(aRow)            
    return sub_matrix    

def det_2_by_2(matrix):
    return matrix[0][0]*matrix[1][1]-matrix[1][0]*matrix[0][1]

def recursion(matrix,somme=None,prod=1):
    if(somme==None):
        somme=[]
    if(len(matrix)==1):
        somme.append(matrix[0][0])
    elif(len(matrix)==2):
        somme.append(det_2_by_2(matrix)*prod)
    else:
        for index, elmt in enumerate(matrix[0]):
            transposee = [list(a) for a in zip(*matrix[1:])]
            del transposee[index]
            mineur = [list(a) for a in zip(*transposee)]
            somme = recursion(mineur,somme,prod*matrix[0][index]*(-1)**(index+2))
    return somme

def M2Maddition(A, B):    
    if(isSquare(A) and isSquare(B)):
        result = zeros_matrix(len(A), len(A))
        for i in range(len(A)):
            for j in range(len(A)):
                result[i][j] = A[i][j] + B[i][j]        
        return result
    else:
        print ("can only be square matrices")   

def norm(m):
    sum = 0
    for i in range(len(m)):
        for j in range(len(m)):
            sum = sum +m[i][j]**2
    return math.sqrt(sum)

def matrix_matrix_subtraction(A, B):    
    if(isSquare(A) and isSquare(B)):
        result = zeros_matrix(len(A), len(A))
        for i in range(len(A)):
            for j in range(len(A)):
                result[i][j] = A[i][j] - B[i][j]        
        return result
    else:
        print ("can only be square matrices") 


def MbyM(m1, m2):    
    if(len(m1[0])==len(m2)):
        result = zeros_matrix(len(m1), len(m2))        
        for i in range(len(m1)):            
            for j in range(len(m2[0])):
                for k in range(len(m1)):
                    result[i][j] +=  m1[i][k] * m2[k][j]
    return result 



def SiedelThis(A, b, xNOT, loops = 20):
    result = []
    result.append(xNOT) # append xNOT
    L = ML(A)
    U = MU(A)
    D = diagonial_matrix(A)

    D_min_L = matrix_matrix_subtraction(D, L)        
    TEMP = M_inverse(D_min_L)
    C = MbyM(TEMP, U)
    if (norm(C)< 1):
        for i in range(loops):
            tm1 = MatrixMult(TEMP, b)
            tm2 = MatrixMult(C, result[i]) 
            x_k1 = V2Vadd(tm1,tm2)
            result.append(x_k1)    
        for index in range(len(result[i])):
            print("X%s = %s"%(index, result[i][index]))
    else:
        print ("does not convergence")

def JacobiMethods(A, b, xNOT, loops = 20):
    result = []
    result.append(xNOT) # append xNOT
    L = ML(A)
    U = MU(A)
    D = diagonial_matrix(A)

    D_inverse = M_inverse(D)        
    TEMP = M2Maddition(L, U)
    C = MbyM(D_inverse, TEMP)
    if (norm(C)< 1):
        for i in range(loops):
            tm1 = MatrixMult(D_inverse, b)
            tm2 = MatrixMult(C, result[i]) 
            x_k1 = V2Vadd(tm1,tm2)
            result.append(x_k1)    
        for index in range(len(result[i])):
            print("X%s = %s"%(index, result[i][index]))
    else:
        print ("does not convergence")

####################################### END OF GAUSS JACOBI ############################################## 

                                                              
                                                              
####################################### NEWTON RAPHSON METHOD ##############################################                                                              
def NewtonRaphson(coeff,dfdx,x,convergence,degree):
    k=0 
    if dfdx(x,coeff)==0:
        print("Halts!!! dfdx = 0")
        return
    print("X%d =%.5f , f(X%d)=%.5f , f'(X%d)=%.5f" % (k,x,k,f(x,degree,coeff),k,dfdx(x,coeff))) 
    while(abs(f(x,degree,coeff))>convergence): 
        k +=1 
        x=x-float(f(x,degree,coeff))/dfdx(x,coeff) 
        print("X%d =%.5f , f(X%d)=%.5f , f'(X%d)=%.5f" % (k,x,k,f(x,degree,coeff),k,dfdx(x,coeff))) 
         
    print ("The value of X=%.5f" %(x)) 
def f(X,degree,coeff): 
    Equation=0 
    indexofCoeff=0 
    for e in coeff: 
        Equation=Equation + e*pow(X,(degree-indexofCoeff)) 
        indexofCoeff +=1 
    return Equation 

def dfdx(X,coeff): 
    Equation2=0 
    indexofCoeff=0 
    for q in range(len(coeff),1,-1): 
        Equation2= Equation2 + (q-1)*coeff[indexofCoeff]*pow(X,q-2) 
        indexofCoeff +=1 
    return Equation2 
    
    
    
    
####################################### END OF NEWTON RAPHSON METHOD ############################################## 


####################################### bisection METHOD ##############################################






#defining a function named Bisection
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 02:22:47 2017

@author: FALALU
"""

#Bisection method
#defining a function named Bisection
def Bisection(F,P1,P2,convergency,degree,coeff):
    aggregate = []
    Intermediate = []
    Xm_Next = []
#    cXmHolder = []
    j = 0
    withHeld1 = float(F(P1,degree,coeff))
    withHeld2 = float(F(P2,degree,coeff))
    
    #checking for sign changes
    if(withHeld1 * withHeld2)>0: 
       print("No solution eist for the inputed equation ")
    Xm = (P1 + P2) / 2.0
    Xu = F(Xm,degree,coeff)
    
    #checking condition for the convergence
    while(abs(Xu) > convergency):
        Xm_Next.append(Xu)
        if((withHeld1 * Xu)>0):
            withHeld1 = Xu
            P1 = Xm
            Intermediate.append(Xm)
        else:
            P2 = Xm
            aggregate.append(Xm)
        Xm=float(P1+P2)/2
        Xu=F(Xm,degree,coeff) 
    for e in  aggregate:
        print("The Values of F(X) within these points is: X%d=%.5f" % (j,e))
        j += 1
    print("The final root of the function is:  %.5f" %(e))
       
def F(X,degree,coeff):
    function = 0
    for e in range(degree + 1):
        function = function + coeff[e] * X**(degree-e)
    return function


def poly():
    guess = []
    f = []
    coeff = []
    degree = int(input("\nEnter the degree of the equation: "))
    for i in range(degree+1):
        k = 1 + i
        coeff.append(float(input("\nEnter the value of NO. %d coefficient of X: " % (k))))
    min = int(input("\nEnter Minimum range: "))
    max = int(input("\nEnter Maximum range: "))    
    
    u = 0
    for j in range(min,max):
        sum = 0
        for k in range(degree + 1):
            sum = sum + coeff[k] * j**(degree - k)
        f.append(sum)

        if(u>0):
            sign = (f[u-1]*f[u])
            if(sign <= 0):
                guess.append(j)
        u = u + 1
    
    print("\nThe guess points on X-axis is/are", guess)
    
    for x in range(len(guess)):
        X1 = int(input("\nChoose the two point values of X from this:  "+ str(guess)+":  "))
        X2 = int(input("\nChoose the two point values of X from this:  "+ str(guess)+":  ")) 
        convergency = float(input("\nEnter the maximum acceptable error: ")) 
    
        solution = Bisection(F,X1,X2,convergency,degree,coeff)
        print(solution)

####################################### bisection METHOD ############################################## 

####################################### SIMULATANEOUS EQUATION SOLVER ############################################## 
#This is a temporary script file.

####################################### END OF SIMULATANEOUS EQUATION SOLVER ############################################## 
        
########################################SIULATANEOUS DETAILS GRABBER #######################################        
#simultaneous eqautions

def solve(matrix, mul=1):
    width = len(matrix)
    if width == 1:
        return mul * matrix[0][0]
    else:
        sign = -1
        total = 0
        for i in range(width):
            m = []
            for j in range(1, width):
                buff = []
                for k in range(width):
                    if k != i:
                        buff.append(matrix[j][k])
                m.append(buff)
            sign *= -1
            total += mul * solve(m, sign * matrix[0][i])
        return total



#end of my function
def Simultaneous():    
    coeff=[]
    n = int(input(" How many equations ? "))
    for i in range(n):
        temp=[]
        print("getting coefficients of equation",i+1)
        for j in range(n):
            r = float(input("coefficient %d of equation %d : "%((j+1),(i+1))))
            temp.append(r)
        coeff.append(temp)
        temp=[]
    re=[]
    for i2 in range(n):
        r2 = float(input("constant value of equation %d : "%(i2+1)))
        re.append(r2)
    det = solve(coeff,1)
#    print(coeff, " is a matrix and the determinant is :",det)
    if det == 0:
        print("no solution")
    else:
        while 1==1:
            print("<<<<<<<<<-------- SIMULTANEOUS EQUATIONS -------->>>>>>>>>>>")
            print("CHOOSE A METHOD TO FIND THE SOLUTION(S)")
            print("Press [1] for CRAMMERS RULE ")
            print("Press [2] for learn about  Gauss Seidel Method")
            print("Press [3] for learn about  Gauss Jacobi Method")
            print("Press [0] to the MAIN MENU")
            CHOICE = int(input("PLEASE ENTER YOUR CHOICE ? :"))
            if CHOICE == 0:
                break
            elif CHOICE == 1:
                Crammers(coeff,re,det,n)
            elif CHOICE == 2:
                SiedelThis(coeff, re, [0,0,0,0], loops = 20)
            elif CHOICE == 3:
                JacobiMethods(coeff, re, [0,0,0,0], loops = 20)
            else:
                print("\nINVALID CHOICE, [0...2] ONLY !!!")
     
        
def Crammers(matrix,constants,det,n):        
        for st in range(n):
            import copy
            matrixN = copy.deepcopy(matrix) # nan wurin fa akwai wahala ashe-ashe baya copy
#            print(matrix," is copied to",matrixN)
            for st2 in range(n):
                matrixN[st2][st]=constants[st2]
#            print(matrix," undergoes change to :",matrixN)
            print("X",st,"=",solve(matrixN,1)/det)
            
#print(len(x))
#print(eq,"\n",const)
    
    

                                                                                                            
####################################### END OF SIMULATANEOUS EQUATION ############################################## 
def LU():
    row = int(input("Enter the number of rows "))
    col = int(input("Enter the no of culumns plus that of b inclusive\ni.e is not a squared matrix(e.g 4colums for 3 rows): "))
    A = []
    temp=[]
    for i in range(row):
        for j in range(col):
            temp.append(float(input("Coefficient eqtn No.%d %d "%((i+1),(j+1))))   )
        A.append(temp)
        temp=[]
    print(A)
    n = len(A) # Give us total of lines

    #  Extract the b vector
    b = [0 for i in range(n)]
    for i in range(0,n):
        b[i]=A[i][n]
    print("Value of column vector b ")
    print (b)

    #  Fill L matrix and its diagonal with 1
    L = [[0 for i in range(n)] for i in range(n)]
    for i in range(0,n):
        L[i][i] = 1
    print("L after transformation ")
    print(L)


    # Fill U matrix
    U = [[0 for i in range(0,n)] for i in range(n)]
    for i in range(0,n):
        for j in range(0,n):
            U[i][j] = A[i][j]
    print("U before transformation ")
    print(U)

    n = len(U)

    # Find both U and L matrices
    for i in range(0,n): # for i in [0,1,2,..,n]
        #  Find the maximun value in a column in order to change lines
        maxElem = abs(U[i][i])
        maxRow = i
        for k in range(i+1, n): # Interacting over the next line
            if(abs(U[k][i]) > maxElem):
                maxElem = abs(U[k][i]) # Next line on the diagonal
                maxRow = k

        # Swap the rows pivoting the maxRow, i is the current row
        for k in range(i, n): # Interacting column by column
            tmp=U[maxRow][k]
            U[maxRow][k]=U[i][k]
            U[i][k]=tmp

        #  Subtract lines
        for k in range(i+1,n):
            c = -U[k][i]/float(U[i][i])
            L[k][i] = c # (4.4) Store the multiplier
            for j in range(i, n):
                U[k][j] += c*U[i][j] # Multiply with the pivot line and subtract


        #  Make the rows bellow this one zero in the current column
        for k in range(i+1, n):
            U[k][i]=0
    print("L after transformation ")
    print (L)
    print("U after transformation ")
    print (U)
    n = len(L)

    # Perform substitutioan Ly=b
    y = [0 for i in range(n)]
    for i in range(0,n,1):
        y[i] = b[i]/float(L[i][i])
        for k in range(0,i,1):
            y[i] -= y[k]*L[i][k]

    print("y values")
    print(y)

    n = len(U)
    print (n)

    # Perform substitution Ux=y
    x = [0 for i in range(n)] ##-1
    for i in range(n-1,0,-1):
    #for i in range (n, -1 , -1):
        print(i)
        x[i] = y[i]/float(U[i][i])
        for k in range (i,0,-1):
            x[i] -= x[k]*U[i][k]
    print (x)
                                                                                                                
###################################### POLYNOMIAL EQUATIONS #############################################
#definition of a polynomial function

def Poly():
    min = int(input("Please Enter the MIN VALUE :"))
    max = int(input("Please Enter the MAX VALUE :"))

    coeff=[];f=[];guess=[];u=0
    degree = int(input("Enter the DEGREE : "))
    for i in range(degree+1):
        #getting the coefficients
        coeff.append(float(input("Enter coeffecient # %d : "%(i+1))))
    
    for j in range(min,max+1):
        #for k in coeff
        sum=0
        for k in range(degree+1):
            sum = sum + coeff[k]*j**(degree-k)
        f.append(sum)
        if(u>1):
            if(f[u-1]*f[u])<=0:
                guess.append(j)
        u=u+1
        
        print("f(%f) = %f "%(j,sum))
        
#    print(guess)
    if len(guess)>0:    
        print("GUESS(ES) are: ")
        for g in guess:
            print(g)
            
        while 1==1:
            print("<<<<<<<<<-------- POLYNOMIAL EQUATIONS -------->>>>>>>>>>>")
            print("CHOOSE A METHOD TO FIND THE SOLUTION(S)")
            print("Press [1] for Newton Raphson Method")
            print("Press [2] for Bisection Method")
            print("Press [0] to the MAIN MENU")
            CHOICE = int(input("PLEASE ENTER YOUR CHOICE ? :"))
            if CHOICE == 0:
                break
            elif CHOICE == 1:
                print("\n<<<-- TESTING WITH ALL THE INITIAL GUESS(ES)-->>>")
                precision = int(input("ENTER THE PRECISIONS i.e., how many decimal points? :"))                
#                for g in guess:
                NewtonRaphson(coeff,dfdx,guess[0],1/(10**precision),degree) 
            elif CHOICE == 2:
                poly()
            else:
                print("\nINVALID CHOICE, [0...2] ONLY !!!")
    else:
        print("THERE IS NO SOLUTION(S) THE GIVEN EQUATION !")
        


###################################### END POLYNOMIAL EQUATIONS #############################################

###################################### LINEAR EQUATIONS #############################################

def Linear():
    while 1==1:
        print("<<<<<<<<<-------- NON LINEAR EQUATIONS -------->>>>>>>>>>>")
        print("Press [1] for SIMUTANEOUS equation")
        print("Press [0] to the MAIN MENU")
        CHOICE = int(input("PLEASE ENTER YOUR CHOICE ? :"))
        if CHOICE == 0:
            break
        elif CHOICE == 1:
                Simultaneous()
        else:
            print("\nINVALID CHOICE, [0...2] ONLY !!!")
    
###################################### LINEAR EQUATIONS #############################################







###################################### NON-LINEAR EQUATIONS #############################################

def NonLinear():
    while 1==1:
        print("<<<<<<<<<-------- NON LINEAR EQUATIONS -------->>>>>>>>>>>")
        print("Press [1] for Quadratic equation")
        print("Press [2] for Polynimial equation\n\t -->Bisection\n\t-->Newton Raphson")
        print("Press [0] to the MAIN MENU")
        CHOICE = int(input("PLEASE ENTER YOUR CHOICE ? :"))
        if CHOICE == 0:
            break
        elif CHOICE == 1:
                Quadratic()
        elif CHOICE == 2:
                Poly()
        else:
            print("\nINVALID CHOICE, [0...2] ONLY !!!")
    
###################################### NON-LINEAR EQUATIONS #############################################







############################  MAIN ENTRY POINT ###############################
#complete menu
def startMenu():
	while 1==1:
		print("<<<<-----------MAIN MENU----------->>>>")
		print("Press [1] for Perfect EGG solutions :")
		print("Press [2] for Linear Equations :")
		print("\t->Crammers Rule ")
		print("\t->Gauss-Seidel Iteration ")
		print("\t->Gauss-Jacobi Iteration ")
		print("Press [3] for Non-Linear Equations :")
		print("\t->Quadratic Equations ")
		print("\t->Polynomial ")
		print("\t\t->Bisection Method")
		print("\t\t->NewtonRapson Method")        
		print("Press [4] for LU Decomposition :")
		print("Press [5] for AUST CGPA Computations :")
		print("Press [0] TO QUIT :")
		CHOICE=int(input("PLEASE ENTER YOUR CHOICE :"))
		if CHOICE == 0:
			break
		elif CHOICE == 1:
			print("\n----->>>>>>>>>>>BOIL PERFECT EGG MODULE<<<<<<<<<<<-----")
			PerfectEgg()
		elif CHOICE == 2:
			print("\n----->>>>>>>>>>>LINEAR EQUATIONS MODULE<<<<<<<<<<<-----")
			Linear()        
		elif CHOICE == 3:
			print("\n----->>>>>>>>>>>NON-LINEAR EQUATIONS MODULE<<<<<<<<<<<-----")
			NonLinear()        
		elif CHOICE == 4:
			LU()
		elif CHOICE == 5:
			print("\n----->>>>>>>>>>>AUST GRADING POINT MODULE<<<<<<<<<<<-----")
			GetGrade()  #ENTRY POINTS 
		else:
			print("\n",CHOICE," is INVALID CHOICE!!!\nPLEASE ENTER 0...9 Only\n")
	print("\n<<<<<<<<<<<<----------GOOD BYE---------->>>>>>>>>>>>\n")
startMenu()
#####################################END OF MAIN ENTRY POINT#########################################
            

