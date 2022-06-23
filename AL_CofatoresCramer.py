# https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
# https://matplotlib.org/3.5.0/tutorials/introductory/pyplot.html

import numpy as np
import math
import cmath
import matplotlib.pyplot as plt

prt=False
pivot=False
A = np.array([[2,1,-3],[-1,3,2],[3,1,-3]])
b = np.array([-1,12,0])

def det_cofat(A):
    N = len(A)
    if N==1:
        return A
    cof = np.zeros(N)
    for col in range(N):
        D = np.delete(A, 0,   0)
        D = np.delete(D, col, 1)
        cof[col] = ((-1)**(col))*det_cofat(D);
    d = A[0,:]@cof
    return d

def det_cofat2(A):
    N = len(A)
    if N==2:
        return A[1,1]*A[2,2]-A[1,2]*A[2,1]
    cof = np.zeros(N)
    for col in range(N):
        D = np.delete(A, 0,   0)
        D = np.delete(D, col, 1)
        cof[col] = ((-1)**(col))*det_cofat(D);
    d = A[0,:]@cof
    return d

def inv_cofat(A,prt):
    N = len(A)
    Cof = np.zeros((N,N))
    for col in range(N):
        for lin in range(N):
            Dij = np.delete(A, lin,   0)
            Dij = np.delete(Dij, col, 1)
            Cof[lin,col] = (-1)**(lin+col)*det_cofat(Dij);
    Adj=Cof.transpose();
    Ai=(1/det_cofat(A))*Adj;
    return Ai

def ResolucaoDireta(A,b):
    x=inv_cofat(A,prt) @ b  
    return x

def Cramer(A,b,prt):
    N = len(A)
    C=np.c_[A, b]
    D=det_cofat(A)
    x = np.zeros(N)
    if (prt): 
        print("Matriz Aumentada [C=A|b] det(A)=%f" % (D))
        print(C)
    for k in range(0,N):
        Ak=A.copy()
        Ak[:,k]=b #substitui coluna k por vetor b
        Dk=det_cofat(Ak)
        x[k]=Dk/D
        if (prt):
            print("Matriz A%d   det(A%1d)=%f" % (k,k,Dk))
            print(Ak)
            print("x(%d)=%f/%f=%f\n" % (k,Dk,D,x[k]))
    return x
