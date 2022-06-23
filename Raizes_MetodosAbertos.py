# https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
# https://matplotlib.org/3.5.0/tutorials/introductory/pyplot.html

import numpy as np
import matplotlib.pyplot as plt
from   Raizes_MetodosIntervalares import *
f = lambda x: 4.9 + (5.1*x-22.3)/(3.1*x-4.2)
f = lambda x: 4*np.exp(x/3)-20*np.exp(-x/5)*np.sin(x)
f = lambda x: x**2-8*x+6
f = lambda x: x**5-6*x**4+10*x*3-10*x**2+5*x-2
#print(np.array([x,f(x)]).transpose())

def  NewtonRaphson(f,df,x0,tol,prt):
    erro = 1;
    if (prt):
        if (x0.imag!=0):
            print ( '%i\t%.10f+%10fi\t%.1e' % (1,x0.real,x0.imag,erro) )
        else:
            print ( '%i\t%.10f\t%.1e' % (1,x0,erro) )    
    for k in range(2,500): 
       f0=f(x0)
       df0=df(x0)
       if (f0==0):
           x1=x0; 
       else:
           x1 = x0 - f0/df0 
       if (x1!=0):
           erro =abs((x1-x0)/x1)
       if (prt): 
            if (x1.imag!=0):
                print ( '%i\t%.10f+%10fi\t%.1e' % (k,x1.real,x1.imag,erro) )
            else:
                print ( '%i\t%.10f\t%.1e' % (k,x1,erro) )   
       if  (erro<tol):
           break
       x0=x1
    return x1


def  NewtonRaphson_p(f,df,x0,tol,p,prt):
    erro = 1;
    if (prt):
        if (x0.imag!=0):
            print ( '%i\t%.10f+%10fi\t%.1e' % (1,x0.real,x0.imag,erro) )
        else:
            print ( '%i\t%.10f\t%.1e' % (1,x0,erro) ) 
    for k in range(2,500): 
       f0=f(x0)
       df0=df(x0)
       if (f0==0):
           x1=x0; 
       else:
           x1 = x0 - p*f0/df0 
       if (x1!=0):
           erro =abs((x1-x0)/x1)
       if (prt): 
            if (x0.imag!=0):
                print ( '%d\t%.10f+%10fi\t%.1e' % (k,x1.real,x1.imag,erro) )
            else:
                print ( '%d\t%.10f\t%.1e' % (k,x1,erro) )     
       if  (erro<tol):
           break
       x0=x1
    return x1

def  NewtonRaphson_pol(ps,x0,prt):
    erro = 1;
    if (prt):
        if (x0.imag!=0):
            print ( '%i\t%.10f+%10fi\t%.1e' % (1,x0.real,x0.imag,erro) )
        else:
            print ( '%i\t%.10f\t%.1e' % (1,x0,erro) )    
    for k in range(2,500): 
       f0=ps(x0)
       ps1 = ps.deriv(1)
       df0=ps1(x0)
       if (f0==0):
           x1=x0; 
       else:
           x1 = x0 - f0/df0 
       if (x1!=0):
           erro =abs((x1-x0)/x1)
       if (prt): 
            if (x1.imag!=0):
                print ( '%i\t%.10f+%10fi\t%.1e' % (k,x1.real,x1.imag,erro) )
            else:
                print ( '%i\t%.10f\t%.1e' % (k,x1,erro) )   
       if  (erro<1e-16):
           print(x1)
           break
       x0=x1
    return x1

def Secante(f,x1,x2,tol,prt):
    erro = 1;
    if (prt):
        print ( '%i\t%.10f\t%.1e' % (1,x1,erro) )
        print ( '%i\t%.10f\t%.1e' % (1,x2,erro) ) 
    for k in range(2,500): 
        x0=x1;
        x1=x2;    
        f0=f(x0)
        f1=f(x1)
        df_x1 = (f1-f0)/(x1-x0);
        x2 = x1 - f(x1)/df_x1;
        if(x2!=0):
            erro =abs((x2-x1)/x2)
        if (prt):
            print ( '%i\t%.10f\t%.1e' % (k,x2,erro) ) 
        if  (erro<tol):
            break
    return x1


def SecanteModificado(f,x0,tol,prt):
    erro = 1;
    dx=1e-8;
    if (prt):
        if (x0.imag!=0):
            print ( '%i\t%.10f+%10fi\t%.1e' % (1,x0.real,x0.imag,erro) )
        else:
            print ( '%i\t%.10f\t%.1e' % (1,x0,erro) )  
    for k in range(2,500): 
       f0=f(x0)
       f1=f(x0+dx)
       df0=(f1-f0) /dx  # Secante Modificado
       x1 = x0 - f0/df0
       if(x1!=0):
           erro =abs((x1-x0)/x1) 
       if (prt): 
            if (x0.imag!=0):
                print ( '%i\t%.10f+%10fi\t%.1e' % (k,x1.real,x1.imag,erro) )
            else:
                print ( '%i\t%.10f\t%.1e' % (k,x1,erro) )  
       if  (erro<tol):
           break
       x0=x1
    return x1

def PontoFixo(f,x0,tol,prt):
    g= lambda x: x + f(x)
    erro = 1;
    if (prt): 
        print ( '%i\t%.10f\t%.1e' % (1,x0,erro)) 
    for k in range(2,500): 
        x1 = g(x0);
        if(x1!=0):
            erro =abs((x1-x0)/x1)
        if (prt): 
            print ( '%i\t%.10f\t%.1e' % (k,x1,erro))
        if  (erro<tol):
            break
        x0=x1
    return x1

#f = lambda x: 4*np.exp(x/3)-20*np.exp(-x/5)*np.sin(x)
#df = lambda x: 4/3*np.exp(x/3)+4*np.exp(-x/5)*np.sin(x)-20*np.exp(-x/5)*np.cos(x)
#x0=2
f = lambda x: x**3-5*x**2+33*x-29
df = lambda x: 3*x**2-10*x+33
x0=10+10j
tol = 0.0001
prt = True
#f = lambda x: np.exp(-x/5)*np.cos(x-1)
