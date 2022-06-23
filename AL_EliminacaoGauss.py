# https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
# https://matplotlib.org/3.5.0/tutorials/introductory/pyplot.html

import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
from tabulate import tabulate

prt=False
pivot=True
A = np.array([[2,1,-3],[-1,3,2],[3,1,-3]])
b = np.array([-1,12,0])

A = np.array([[22,1,-3,6],[-1,3,2,-2],[-4,3,1,-3],[-8,3,2,1]])
b = np.array([-1,12,3,6])


def EliminarLinha(lin,p,C,prt):
  [N,M]=C.shape
  if (prt):
    print("(L%d)=(L%d)-(%f)/(%f)*L(%d)" % (lin,lin,C[lin,p],C[p,p],p))
    if (lin<N):
        print("\n")
  m=C[lin,p]/C[p,p];
  C[lin,p:M+1]=C[lin,p:M+1]-m*C[p,p:M+1];
  return (C,m)


def PivotarColuna(p,l_ini,l_fim,C,prt):
     [N,M]=C.shape
     linhas = np.arange(l_ini,l_fim)
     #max_lin_p=np.amax(abs(C[linhas,p]));  
     dist=np.argmax(abs(C[linhas,p]))
     if(dist!=0):  
         C[[ p, (dist+p)],:] = C[[ (dist+p) , p ],:] #troca linhas
         if(prt):
             print("Trocando linhas %d e %d" % (p,dist+p-1))
             pdtabulate=lambda C:tabulate(C,headers='keys')
             print(pdtabulate(C)) 
     return (C,dist)

def SubstituicaoProgressiva(A,b,prt):
  [N,M]=A.shape
  y = np.zeros(N)
  y[0]=b[0];
  if (prt):
      print("Substituição Progressiva\ny(%d)=%f\n" % (1,y(1)))
  for lin in range(1,N):
        print(b[lin])
        y[lin]= b[lin]-A[lin,0:lin] @ y[0:lin];
        print(y[lin])
        if (prt):
            print("y(%d)=%f\n" % (lin,y(lin)))
  return y

def SubstituicaoRegressiva(A,b,prt):
  [N,M]=A.shape
  x = np.zeros(N)
  x[N-1]=b[N-1]/A[N-1,N-1];
  if (prt):
      print("Substituição regressiva\nx(%d)=%f\n" % (N-1,x[N-1]))
  for lin in range(N-2,-1,-1):    
     x[lin]=(b[lin]-A[lin,lin+1:N] @ x[lin+1:N])/A[lin,lin];
     if (prt):
         print("x(%d)=%f\n" % (lin,x[lin]))
  return x


def EliminacaoGauss(A,b,prt,pivot): # com pivotamento
   A= A.astype(float)
   b=b.astype(float)
   C=np.c_[A, b]
   [N,M]=C.shape;
   x = np.zeros(N)
   if(prt):
       print("Matriz Aumentada [C=A|b]")
       print(C)
   for p in range(N-1):
       if (pivot):
         (C,dist)=PivotarColuna(p,p,N,C,prt)
       elif (C[p,p]== 0):
          break
       if (prt):
            print("Eliminando coluna %d com Pivô %f\n" % (p,C[p,p]) )
       for lin in range(p+1,N): #eliminação progressiva
           (C,m)=EliminarLinha(lin,p,C,prt)
       if (prt):
           print(C)
   if (C[p,p]!=0): 
         x=SubstituicaoRegressiva(C[:,0:N],C[:,N],prt)
   else:
         print("Não há solução única pois matriz A é singular")
         x[1:N]= np.inf
   return x

def det_gauss(A,prt,pivot): # com pivotamento
   A= A.astype(float)
   [N,M]=A.shape;
   if(prt):
       print("Matriz A")
       print(A)
   npivot=0;    
   for p in range(N-1):
       if (pivot):
         (A,dist)=PivotarColuna(p,p,N,A,prt)
         if(dist!=1):
               npivot=npivot+1
       elif (A[p,p]== 0):
          break
       if (prt):
            print("Eliminando coluna %d com Pivô %f\n" % (p,A[p,p]) )
       for lin in range(p+1,N): #eliminação progressiva
           (A,m)=EliminarLinha(lin,p,A,prt)
       if (prt):
           print(A)
   d=np.prod(np.diag(A)) * (-1)**npivot
   return d


def EliminacaoGaussJordan(A,b,prt,pivot): # com pivotamento
   A= A.astype(float)
   b=b.astype(float)
   C=np.c_[A, b]
   [N,M]=C.shape;
   x = np.zeros(N)
   if(prt):
       print("Matriz Aumentada [C=A|b]")
       print(C)
   for p in range(N):
       if (pivot):
         (C,dist)=PivotarColuna(p,p,N,C,prt)
       elif (C[p,p]== 0):
          break
       C[p,:]=C[p,:]/C[p,p]; 
       if (prt):
            print("Eliminando coluna %d com Pivô %f\n" % (p,C[p,p]) )
       linhas = np.concatenate((np.arange(0,p),np.arange(p+1,N)))
       for lin in linhas: #eliminação progressiva
           (C,m)=EliminarLinha(lin,p,C,prt)
       if (prt):
           print(C)  
   if (C[p,p]!=0): 
         x=C[:,N]
   else:
         print("Não há solução única pois matriz A é singular")
         x[1:N]= np.inf
   return x


def inv_gauss(A,prt,pivot):
   A= A.astype(float)
   [N,M]=A.shape;
   C=np.c_[A, np.eye(N)]
   [N,M]=C.shape;
   Ai = np.zeros((N,N))
   if(prt):
       print("Matriz Aumentada [C=A|Ai]")
       print(C)
   for p in range(N):
       if (pivot):
         (C,dist)=PivotarColuna(p,p,N,C,prt)
       elif (C[p,p]== 0):
          break
       C[p,:]=C[p,:]/C[p,p]; 
       if (prt):
            print("Eliminando coluna %d com Pivô %f\n" % (p,C[p,p]) )
       linhas = np.concatenate((np.arange(0,p),np.arange(p+1,N)))
       for lin in linhas: #eliminação progressiva
           (C,m)=EliminarLinha(lin,p,C,prt)
       if (prt):
           print(C)
   if (C[p,p]!=0): 
         Ai=C[0:N,N:2*N]
   else:
         print("Não há solução única pois matriz A é singular")
         Ai= []
   return Ai

def FatoracaoLU(A,prt,pivot): # com pivotamento
   A= A.astype(float)
   [N,M]=A.shape;
   L=np.eye(N)
   P=L.copy()
   U=A.copy()
   if(prt):
       print("Matriz Aumentada [C=A|b]")
       print(np.c_[L, U])
   for p in range(N-1):
       if (pivot):
         (U,dist)=PivotarColuna(p,p,N,U,prt)
       if(dist!=1):
          P[[ p, (dist+p)],:] = P[[ (dist+p) , p ],:]
          for col in range(p): # troca linhas abaixo da diagonal
            if (col<p):
               L[[p,(dist+p)],col] = L[[(dist+p),p],col] 
          if (prt):
            print("Trocando linhas %d e %d de U\n" % (p,p+dist)) 
            print("Trocando linhas %d e %d de L abaixo da diagonal" % (p,p+dist))
            print(np.c_[L, U])
       if (U[p,p]==0):
           break;
       for lin in range(p+1,N): #eliminação progressiva
           (U,m)=EliminarLinha(lin,p,U,prt)
           L[lin,p] = m
       if (prt):
           print(np.c_[L, U])
   return (L,U,P)

def SolucaoLU(L,U,P,b,prt):
    bp=P @ b
    y=SubstituicaoProgressiva(L,bp,False)
    x=SubstituicaoRegressiva(U,y,False)
    if (prt):
        print("[b bp]")
        print(np.c_[b, bp])
        print("matriz aumentada  [L  bp y]")
        print(np.c_[L, bp, y])
        print("matriz aumentada  [U y x]")
        print(np.c_[U, y, x])
    return x

dl=np.array([7,9,4,1,2,34,32,65,23,12,12,5,9,3]);

dp=np.array([1,3,5,3,4,6,7,5,3,5,6,8,10,7,14]);

du=np.array([3,7,3,2,1,5,3,9,8,24,12,6,7,8]);

r=np.array([1,1.98,4.99,2.98,6.97,4.01,4.98,2.25,3.37,6.94,4.92,3.04,5.97,7.91,1.89]);

def tridiagonal(dl,dp,du,r):
    dl= dl.astype(float)
    dp= dp.astype(float)
    du= du.astype(float)
    r= r.astype(float)
    N=len(r)
    y = np.zeros(N)
    for k in range(1,N): #eliminação progressiva
        m=dl[k-1]/dp[k-1];
        dp[k]=dp[k]-m*du[k-1];
        r[k]=r[k]-m*r[k-1];
    y[N-1]=r[N-1]/dp[N-1]; #substituição regressiva
    for k in range(N-2,-1,-1):
        y[k]=(r[k]-du[k]*y[k+1])/dp[k];
    return y


#np.concatenate((np.arange(3,6),np.arange(1,2)))
