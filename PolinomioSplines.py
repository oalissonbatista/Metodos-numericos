import numpy as np
import matplotlib.pyplot as plt
import io
from scipy import linalg
from   AL_EliminacaoGauss import *
from tabulate import tabulate

x=np.array([1,2,4,7,11]);
y=np.array([2,5,7,6,1]);

def  SplineLinear(x,y): 
     N=len(x)-1 # de intervalos 
     b   =   np.zeros([2*N]);
     A   =   np.zeros([2*N,2*N]);
     ps  =   []
     index=0; # x(index),y(index)
     col=0;   # A(lin,col)
     for lin in range (0,2*N,2):  #2N equações: funções passam por (xi,yi))
        A[lin:lin+2,col:col+2]= [[x[index],   1],
                                 [x[index+1], 1]]  
        b[lin:lin+2]=  [y[index],
                        y[index+1]] ;
        index=index+1;
        col=col+2;  
     coef=tridiagonal(np.diag(A,-1),np.diag(A),np.diag(A,1),b)
     for i in range (N):
        u=np.array([coef[2*i],coef[2*i+1]])
        ps.append(np.poly1d(u)) 
     return (ps,A,b)

def  SplineQuadratica(x,y): 
     N=len(x)-1 # de intervalos 
     b   =   np.zeros([3*N]);
     A   =   np.zeros([3*N,3*N]);
     ps  = []
     index=0;  # x(index),y(index)
     col=0;    # A(lin,col)
     A[0,0]=1; # 1a equação a1=0 
     for lin in range (1,2*N,2):   # 2N equações: funções passam por (xi,yi))
        A[lin:lin+2,col:col+3]= [[x[index]**2,   x[index],   1],
                                 [x[index+1]**2, x[index+1], 1]]  
        b[lin:lin+2]=  [y[index],
                        y[index+1]] ;
        index=index+1;
        col=col+3;
     col=0;
     index=1;
     for lin in range (2*N+1,3*N): # N-1 equações-1a derivadas pontos internos
        A[lin,col:col+6]=   [2*x[index], 1, 0, -2*x[index], -1, 0];
        index=index+1;
        col=col+3;
     #coef=np.linalg.solve(A,b)
     coef=EliminacaoGaussJordan(A,b,False,True);
     for i in range (N):
         u=np.array([coef[3*i],coef[3*i+1],coef[3*i+2]])
         ps.append(np.poly1d(u)) 
     return (ps,A,b)

def  SplineCubica(x,y): 
     N=len(x)-1 # de intervalos 
     b   =   np.zeros([4*N]);
     A   =   np.zeros([4*N,4*N]);
     ps  = []
     index=0;  # x(index),y(index)
     col=0;    # A(lin,col)
     for lin in range (0,2*N,2):   # 2N equações: funções passam por (xi,yi))
        A[lin:lin+2,col:col+4]= [[x[index]**3,   x[index]**2,   x[index],   1],
                                 [x[index+1]**3, x[index+1]**2, x[index+1], 1]]  
        b[lin:lin+2]=  [y[index],
                        y[index+1]] ;
        index=index+1;
        col=col+4;
     col=0;
     index=1;
     for lin in range (2*N,3*N-1):   # N-1 equações - 1a derivadas pontos internos
        A[lin,col:col+8]=[3*x[index]**2,2*x[index],1,0,-3*x[index]**2,-2*x[index],-1,0];
        index=index+1;
        col=col+4;
     col=0;
     index=0;  
     A[3*N-1,col:col+4]=   [6*x[index],   2,  0, 0]; #2a derivada ponto externo 
     index=index+1;
     for lin in range (3*N,4*N-1):   # N-1 equações - 2a derivadas pontos internos
        A[lin,col:col+8]=   [6*x[index],   2,  0, 0, -6*x[index], -2, 0, 0];
        index=index+1;
        col=col+4;
     A[4*N-1,col:col+4]=   [6*x[index],   2,  0, 0]; #2a derivada ponto externo
     coef=EliminacaoGaussJordan(A,b,False,True);
     #coef=np.linalg.solve(A,b)
     for i in range (N):
         u=np.array([coef[4*i],coef[4*i+1],coef[4*i+2],coef[4*i+3]])
         ps.append(np.poly1d(u)) 
     return (ps,A,b)
   
def SplineNatural(x,y):
    N=len(x)-1  # de intervalos
    h  =np.diff(x);
    dy =np.diff(y)/h;
    dl=h[1:N-1];
    dp=2.0*(h[0:N-1]+h[1:N]);
    du=h[1:N-1];
    r=3.0*(dy[1:N]-dy[0:N-1]);
    a=y;
    b=np.zeros(N)
    c = np.concatenate([[0.0],tridiagonal(dl,dp,du,r), [0.0]]);
    d=np.zeros(N)
    ps=[]
    for i in range(N):
      b[i] = dy[i] - 1/3*h[i]*(2*c[i]+c[i+1]);
      d[i] = (c[i+1]-c[i])/(3*h[i]);
      u=np.array([d[i],c[i],b[i],a[i]])
      ps.append(np.poly1d(u))
    return (ps,a,b,c,d,dl,dp,du,r)
   
def Plot_Spline(x,y,tipo):
      if (tipo==1):
            (ps,A,b)=SplineLinear(x,y)
      elif (tipo==2):
            (ps,A,b)=SplineQuadratica(x,y)
      else:
            (ps,A,b)=SplineCubica(x,y)
      print("Matriz Aumentada [A|b]")
      (N,N)=A.shape
      C = np.insert(A, N, b, axis=1)
      pdtabulate=lambda C:tabulate(C,headers='keys')
      print(pdtabulate(C)) # coef = inv(A)*b
      K=len(ps)
      for i in range (0,K):
          print("ps(%d)=" % (i))
          print(ps[i])
          xp=np.linspace(x[i],x[i+1],50);
          yp=ps[i](xp)
          plt.plot(xp,yp)
      plt.scatter(x,y)
      plt.grid()
      plt.title("Interpolação Splines");
      plt.show()
      return


def Plot_SplineNatural(x,y):
      (ps,a,b,c,d,dl,dp,du,r)=SplineNatural(x,y)
      K=len(ps)
      print("[x h y dy r]");
      c2=np.concatenate([np.diff(x), [np.nan]])
      c4=np.concatenate([np.diff(y)/np.diff(x), [np.nan]])
      c5=np.concatenate([[np.nan], r, [np.nan] ])
      c6=np.concatenate([[np.nan], dp, [np.nan]])
      C=np.array([x,c2,y,c4,c5,c6]).transpose()
      pdtabulate=lambda C:tabulate(C,headers='keys')
      print(pdtabulate(C))
      print("Sistema Tridiagonal [dl dp du | r = c]");
      C=np.diag(dp)+np.diag(du,1)+np.diag(dl,-1)
      C = np.insert(C, K-1, r, axis=1)
      C = np.insert(C, K, c[1:K], axis=1)
      pdtabulate=lambda C:tabulate(C,headers='keys')
      print(pdtabulate(C))
      print("Matriz de Soluções [a b c d]");
      print(np.array([a[0:K],b,c[0:K],d]).transpose())
      print("Splines Cúbicas");
      for i in range (0,K):
          print("ps(%d)=" % (i))
          print(ps[i])
          xp=np.linspace(x[i],x[i+1],50);
          yp=ps[i](xp-x[i])
          plt.plot(xp,yp)
      plt.scatter(x,y)
      plt.grid()
      plt.title("Interpolação Splines");
      plt.show()
      return
