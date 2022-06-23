# https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
# https://matplotlib.org/3.5.0/tutorials/introductory/pyplot.html

import numpy as np
import matplotlib.pyplot as plt
from   PolinomioInterpolacao import *

def Trapezios_tab(x,y,prt,plot):
   n=len(x)
   h=x[1]-x[0] #espaçamento
   w=np.ones(n)
   w[1:-1]=2
   S = y @ w
   I=(h/2)*S
   if (prt):
      print("[x y w y*w]")
      print(np.array([x,y,w,y*w]).transpose())
      print("h=%f\nsoma(y*w)=%f\nI=(%f)/2*(%f)\n" % (h,S,h,S))
      print("Integal por Trapézios %4d pontos = %f\n" % (n,I))
   if (plot):
      for i in range(0,n-1):
         (pv,A,coef)= PolinomioVandermonde(x[i:i+2],y[i:i+2])
         t = np.linspace(x[i],x[i+1],50);
         plt.plot(t,pv(t))
      plt.scatter(x,y)
      plt.xlabel("Integral por Trapézios " + str(n) + " pontos=" + str(I))
      plt.grid(True)
      plt.show()
   return I

def Simpson_13_tab(x,y,prt,plt1):
   n=len(x)
   h=x[1]-x[0] #espaçamento
   w=np.ones(n)
   w[1:-1:2]=4
   w[2:-1:2]=2
   S = y @ w
   I=(h/3)*S
   if (prt):
      print("[x y w y*w]")
      print(np.array([x,y,w,y*w]).transpose())
      print("h=%f\nsoma(y*w)=%f\nI=(%f)/3*(%f)\n" % (h,S,h,S))
      print("Integal por Simpson 1/3 %4d pontos = %f\n" % (n,I))
   if(plt1):
      for i in range(0,n-2,2):
         (pv,A,coef)= PolinomioVandermonde(x[i:i+3],y[i:i+3])
         t = np.linspace(x[i],x[i+2],50);
         plt.plot(t,pv(t))
      plt.scatter(x,y)
      plt.xlabel("Integral por Simpson 1/3 " + str(n) + " pontos=" + str(I))
      plt.grid(True)
      plt.show()
   return I

def Simpson_38_tab(x,y,prt,plt1):
   n=len(x)
   h=x[1]-x[0] #espaçamento
   w=np.ones(n)
   w[1:-1:3]=3
   w[2:-1:3]=3
   w[3:-1:3]=2
   S = y @ w
   I=(3*h/8)*S
   if (prt):
      print("[x y w y*w]")
      print(np.array([x,y,w,y*w]).transpose())
      print("h=%f\nsoma(y*w)=%f\nI=3*(%f)/8*(%f)\n" % (h,S,h,S))
      print("Integal por Simpson 3/8 %4d pontos = %f\n" % (n,I))
   if (plt1):
      for i in range(0,n-3,3):
         (pv,A,coef)= PolinomioVandermonde(x[i:i+4],y[i:i+4])
         t = np.linspace(x[i],x[i+3],50);
         plt.plot(t,pv(t))
      plt.scatter(x,y)
      plt.xlabel("Integral por Simpson 3/8 " + str(n) + " pontos=" + str(I))
      plt.grid(True)
      plt.show()
   return I

def Trapezios(f,a,b,n,prt,plt1):
   x = np.linspace(a,b,n)
   y = f(x)
   Trapezios_tab(x,y,prt,plt1)

def Simpson_13(f,a,b,n,prt,plt1):
   x = np.linspace(a,b,n)
   y = f(x)
   Simpson_13_tab(x,y,prt,plt1)

def Simpson_38(f,a,b,n,prt,plt1):
   x = np.linspace(a,b,n)
   y = f(x)
   Simpson_38_tab(x,y,prt,plt1)

x = np.arange(0,12.5,1.25)
y = [0,0.68,3.44,2.78,5.65,6.15,7.01,9.65,9.01,11.97]
prt = True
plt1 = True
#Simpson_13_tab(x, y, True)
f = lambda x: x + np.sin(3.0*x)
#Simpson_13(f, 0, 10, 9, True)
