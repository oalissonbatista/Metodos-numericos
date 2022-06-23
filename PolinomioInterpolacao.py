import numpy as np
import matplotlib.pyplot as plt
import io
from scipy import linalg
from   AL_EliminacaoGauss import *
from tabulate import tabulate

x=np.array([3.0,4.5,6.0,9.0])
y=np.array([5.0,7.1,8.7,10.0])
A=np.array([[2,1,-3],[-1,3,2],[3,1,-3]])

x=np.array([1,2,4,7,11]);
y=np.array([2,5,7,6,1]);

def PolinomioVandermonde(x,y):
   N = len(x)
   A = np.zeros((N,N))
   for k in range(N):
      A[:,k] = x**k
   c=linalg.solve(A,y)
   pv = np.poly1d(c[-1::-1])
   return (pv,A,c) 

def Plot_PolinomioVandermonde(x,y):
    (pv,A,coef)=PolinomioVandermonde(x,y)
    print("Matriz Aumentada [A|y]")
    (N,N)=A.shape
    C = np.insert(A, N, y, axis=1)
    print(C)
    print("Coeficientes coef=A^-1 * y")
    print(coef)
    print("Polinômio de Vandermonde\npv(s)=")
    print(pv)
    xp=np.linspace(min(x),max(x),1000)
    yp=pv(xp)
    plt.plot(xp,yp)
    plt.scatter(x,y)
    plt.title("Polinômio de Vandermonde");
    plt.show()
    return
   
def  PolinomioNewton(x,y):
    N=len(x);   #Cálculo dos coeficientes
    D=np.zeros([N,N])
    D[:,0]=y;
    for j in range (1,N):  # D(N,N) - Tabela de diferenças
        for k in range (j,N):
            D[k,j]=(D[k,j-1]-D[k-1,j-1])/(x[k]-x[k-j])   
    coef=np.diag(D)
    return (coef,D)

def  InterNewton(xp,x,coef):
  N=len(coef)  
  yp=coef[N-1];  
  for k in range ((N-2),-1,-1):
       print(k)
       yp=yp*(xp-x[k]) + coef[k];
  return yp


def Plot_PolinomioNewton(x,y):
    (coef,D)=PolinomioNewton(x,y) 
    N=len(x) 
    print("Tabela de diferenças")
    print(D)
    print("Coeficientes coef=diag(D)")
    print(coef)
    print("Polinômio de Newton\npn(s)=")
    for k in range (N,1,-1):
        print("(%.5e)" % (coef[N-k+1]), end = ' ')
        for j in range ((N),k-1,-1):
           print("(x-(%.5f))" % (x[N-j]), end = ' ')  
        if (k>2):
           print(" + ")
    xp=np.linspace(min(x),max(x),1000);
    yp=InterNewton(xp,x,coef);
    plt.plot(xp,yp)
    plt.scatter(x,y)
    plt.title("Polinômio de Newton")
    plt.show()
    return;

def InterLagrange(xp,x,y):
    N=len(x)
    M=len(xp)
    yp=np.zeros(M)
    Li=np.zeros(N)
    for k in range (0,M):
       for i in range (0,N):
           Li[i]=y[i]
           for j in range (0,N):
             if (j!=i):
                Li[i]=Li[i]*(xp[k]-x[j])/(x[i]-x[j]); 
       yp[k]=np.sum(Li)
    return yp

def print_to_string(*args, **kwargs):
    output = io.StringIO()
    print(*args, file=output, **kwargs)
    contents = output.getvalue()
    output.close()
    return contents

def Plot_PolinomioLagrange(x,y):
    N=len(x); 
    print("Polinômio de Lagrange")
    linha_s = "------------------------------------------------------"
    for i in range (0,N):
        num_s = ""
        den_s = ""
        for j in range (0,N):
           if(i!=j):
             num_s=print_to_string("%s(x-(%.4f))" % (num_s,x[j]), end = ' ')
             den_s=print_to_string("%s(%.4f)" % (den_s,x[i]-x[j]), end = ' ')
        print("L[%d,%d](x)=\n\t%s\n(%.4f)%s\n\t%s" % (i,N-1,num_s,y[i],linha_s,den_s))
    print("pl=")
    for i in range (0,N-1):
       print("y(%d)*L[%d,%d](x)+" % (i-1,i,N-1) , end = ' ')
    print("y(%d)*L[%d,%d](x)" % (N-1,N-1,N-1))
    xp=np.linspace(min(x),max(x),1000);
    yp=InterLagrange(xp,x,y);
    plt.plot(xp,yp)
    plt.grid()
    plt.scatter(x,y)
    plt.title("Polinômio de Lagrange");
    plt.show()
    return


