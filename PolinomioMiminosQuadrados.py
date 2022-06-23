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


def PolinomioMQ(x,y,ordem):
    N = len(x)
    M=ordem+1
    A = np.zeros((N,M))
    for k in range(M):
        A[:,k]=x**k
    u=linalg.solve(A.transpose()@A,A.transpose()@y)
    pm = np.poly1d(u[-1::-1])
    return (pm,A,u)

def Plot_PolinomioMQ(x,y,ordem):
    (pm,A,u)=PolinomioMQ(x,y,ordem)
    (N,M)=A.shape
    print("Matriz Aumentada [(A' A)|(A' y))]")
    C = np.insert(A.transpose()@A, M, A.transpose()@y, axis=1)
    print(C)
    print("Coeficientes u=inv(A' A)*(A' y))")
    print(u)
    print("Polinômio Mínimos Quadradros ordem %d\npm=" % (ordem));
    print(pm)
    St = sum( (y-np.mean(y))**2 ); #goodness of fitness
    Sr = sum( (y-A@u)**2 );
    r2= abs(St-Sr)/St
    print("r2=%.1f%%\tr=%.3f\n" % (r2*100,np.sqrt(r2)) );
    xp=np.linspace(min(x),max(x),1000);
    yp=pm(xp)
    plt.plot(xp,yp)
    plt.grid()
    plt.scatter(x,y)
    plt.title("Polinômio de Mínimos Quadrados");
    plt.show()
    return

def   MQL(x,y):
      n=len(x)
      A=np.array([[n,sum(x)],[sum(x),sum(x**2)]])
      b=np.array([sum(y),sum(x*y)])
      u=linalg.solve(A,b)
      a0=u[0]
      a1=u[1]
      return(a0,a1, A, b)

def plot_MQL(x,y):
    (a0,a1, A, b)=MQL(x,y)
    St = sum( (y-np.mean(y))**2 );
    Sr = sum( (y-a0-a1*x)**2 );
    r2= abs(St-Sr)/St
    (N,M)=A.shape
    C = np.insert(A, M, b, axis=1)
    print("Matriz Aumentada [(A)|(b))]")
    print(C)
    print("y=%.3f+%.3f*x" % (a0,a1))
    print("r2=%.1f%%\tr=%.3f" % (r2*100,np.sqrt(r2)));
    xp=np.linspace(min(x),max(x),1000);
    plt.plot(xp,a0+a1*xp)
    plt.grid()
    plt.scatter(x,y)
    plt.title("Mínimos Quadrados Linear");
    plt.show()
    return
