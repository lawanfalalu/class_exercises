# -*- coding: utf-8 -*-
"""
Created on Tue May 30 23:21:57 2017

@author: FALALU
"""

import math
m = []
def matrix_L(m):
    if(isSquare(m)):
        zeros_diagonal = zeros_matrix(len(m), len(m))        
        for i in range(1, len(m), 1):
            for j in range(0, i):
                zeros_diagonal[i][j] = - m[i][j]
        return zeros_diagonal
    else:
        print ("Must be a square matrix")

def matrix_U(m):
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

def vector_vector_add(v1, v2):
    result = zeros_array(len(v1))
    for i in range(len(v1)):
        result[i] = v1[i]+v2[i]
    return result

def matrix_inverse(m):
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
    
def matrix_vector_mult(m, v):
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

def matrix_matrix_addition(A, B):    
    if(isSquare(A) and isSquare(B)):
        result = zeros_matrix(len(A), len(A))
        for i in range(len(A)):
            for j in range(len(A)):
                result[i][j] = A[i][j] + B[i][j]        
        return result
    else:
        print ("Must be square matrices")   

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
        print ("Must be square matrices") 


def matrix_matrix_mult(m1, m2):    
    if(len(m1[0])==len(m2)):
        result = zeros_matrix(len(m1), len(m2))        
        for i in range(len(m1)):            
            for j in range(len(m2[0])):
                for k in range(len(m1)):
                    result[i][j] +=  m1[i][k] * m2[k][j]
    return result 



def L_U_Siedel(A, b, x_0, iterations = 20):
    result = []
    result.append(x_0) # append x_0
    L = matrix_L(A)
    U = matrix_U(A)
    D = diagonial_matrix(A)

    D_min_L = matrix_matrix_subtraction(D, L)        
    TEMP = matrix_inverse(D_min_L)
    C = matrix_matrix_mult(TEMP, U)
    if (norm(C)< 1):
        for i in range(iterations):
            tm1 = matrix_vector_mult(TEMP, b)
            tm2 = matrix_vector_mult(C, result[i]) 
            x_k1 = vector_vector_add(tm1,tm2)
            result.append(x_k1)    
        for index in range(len(result[i])):
            print("X%s = %s"%(index, result[i][index]))
    else:
        print ("Non - convergence")

def L_U_Jacobi(A, b, x_0, iterations = 20):
    result = []
    result.append(x_0) # append x_0
    L = matrix_L(A)
    U = matrix_U(A)
    D = diagonial_matrix(A)

    D_inverse = matrix_inverse(D)        
    TEMP = matrix_matrix_addition(L, U)
    C = matrix_matrix_mult(D_inverse, TEMP)
    if (norm(C)< 1):
        for i in range(iterations):
            tm1 = matrix_vector_mult(D_inverse, b)
            tm2 = matrix_vector_mult(C, result[i]) 
            x_k1 = vector_vector_add(tm1,tm2)
            result.append(x_k1)    
        for index in range(len(result[i])):
            print("X%s = %s"%(index, result[i][index]))
    else:
        print ("Non - convergence")
if __name__ == "__main__":
    L_U_Jacobi([[10,1,2,-1],[3,-25,-2,1],[2,-1,-6,2],[2,-3,-0.5,-8]],[1,0,-4,3],[0,0,0,0],20)
