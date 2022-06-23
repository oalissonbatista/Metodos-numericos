# https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
# https://matplotlib.org/3.5.0/tutorials/introductory/pyplot.html

import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
from Raizes_MetodosAbertos import NewtonRaphson_pol 
from FE_FloatingPoint import SomaFloat

prt = False
ps=np.poly1d([ 1, -5,  6])  # x^2-8s+6
p_div=np.poly1d([1,-3])
#ps=np.poly1d([ 1, -60,  1100, -6000])  # x^3-60x^2+1100s-6000
ps=np.poly1d([1,-10,31,-30])
p_div =np.poly1d([ 1, -5,  6])
ps3=np.poly1d([ 1., 11., 23., 45.])

def clean_array(a):
    index =  np.argwhere(np.abs(a) < 1e-12)
    if(len(index)>0):
        a[index] = 0
    index =  np.argwhere(np.abs(a.imag) < 1e-12)
    if(len(index)>0):
        a[index] = a[index].real
    return a

def clean_raiz(a,tol):
    if (np.abs(a) < tol):
        a = 0
    if (np.abs(a.imag) < tol):
        a = a.real
    return a

def baskara(ps):
    a=np.flip(ps.c).real;
    delta=a[1]**2-4*a[2]*a[0];
    r1=(-a[1] + cmath.sqrt(delta))/(2*a[2])
    r2=(-a[1] - cmath.sqrt(delta))/(2*a[2])
    return(r1,r2,delta)

def BaskaraFormula(ps):
    [a0, a1] = (ps[0]/ps[2], ps[1]/ps[2])
    delta=a1**2-4*a0;
    S1 = cmath.sqrt(delta)/2
    k=np.array([1,2])
    rts= -a1/2 + S1*np.exp(2j*math.pi*k/2)
    index =  np.argwhere(np.abs(rts.imag) < 1e-12)
    rts=clean_array(rts)
    return (rts,delta)

def CardanoFormula(ps):
    [a0, a1, a2] = (ps[0]/ps[3], ps[1]/ps[3], ps[2]/ps[3])
    Q=(3*a1-a2**2)/9
    R=(9*a2*a1-27*a0-2*a2**3)/54
    delta=Q**3+R**2
    S1=(R+cmath.sqrt(delta))**(1/3)
    if (Q==0):
        S2=0;
    else:
        S2=-Q/S1; 
    k=np.array([1,2,3])  
    rts=-a2/3+S1*np.exp(2j*math.pi*k/3)+S2*np.exp(-2j*math.pi*k/3)
    rts=clean_array(rts)
    return rts,delta

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
           break
       x0=x1   
    return x1


def fatorar_p1(ps,p_div):
    a=np.flip(ps.c).real
    c=np.flip(p_div.c).real
    n=len(ps)
    b=np.zeros(n+1)
    b[n]=a[n];
    for k in range (n-1,-1,-1):
       b[k]=a[k]-c[0]*b[k+1];
    resto = b[0]
    ps_out=np.poly1d(np.flip(b[1:n+1]))   
    return (ps_out, resto, b)

def fatorar_p2(ps,p_div):
    a=np.flip(ps.c).real
    c=np.flip(p_div.c).real
    n=len(ps)
    b=np.zeros(n+1)
    b[n]=a[n];
    b[n-1]=a[n-1]-c[1]*b[n];
    for k in range (n-2,-1,-1):
        b[k]=a[k]-b[k+1]*c[1]-b[k+2]*c[0];
    resto = np.poly1d([b[1],b[0]+b[1]*c[1]])
    ps_out=np.poly1d(np.flip(b[2:n+1]))
    return (ps_out, resto, b)

def muller_iterativo(ps,prt):
    n=len(ps)+1 # polinônio ordem n-1
    rts=np.zeros(n)*(1+1j)
    while (n>3):  # repita enquanto ordem do polinômio for maior que 2 
          raiz=NewtonRaphson_pol(ps,0.1+0.1j,False);
          rts[n-1]=clean_raiz(raiz,1e-8)
          if(np.abs(rts[n-1].imag) < 1e-8 ):
               divisor =  np.poly1d([1,-rts[n-1]]);
               if (prt):
                   print(ps)
                   print(divisor)
                   print("{:.5f}\n".format(rts[n-1]))
               (ps, resto, b) = fatorar_p1(ps, divisor)   # ps=pdiv(ps,divisor)
               n=n-1;   # diminua em 1 a ordem do polinômio
          else:  # se a raiz for complexa elimine também o conjugado
              rts[n-2]=np.conj(rts[n-1])
              r=rts[n-1].real
              i=rts[n-1].imag
              divisor = np.poly1d([1,-2*rts[n-1].real,abs(rts[n-1].real)**2])
              if (prt):
                  print(ps)
                  print(divisor)
                  print("{:.5f}".format(rts[n-1]))
                  print("{:.5f}\n".format(rts[n-2]))
              (ps, resto, b) = fatorar_p2(ps, divisor)  # ps=pdiv(ps,divisor)
              n=n-2;   # diminua em 2 a ordem do polinômio
    if(n==3):   # último polinômio de ordem 1 ou 2  
        [rts[n-1],rts[n-2],delta]=baskara(ps);  
        if (prt):
            print(ps)
            print("{:.5f}".format(rts[n-1]))
            print("{:.5f}\n".format(rts[n-2]))
    else: 
        rts[n-1]=-ps[0]/ps[1]
        if (prt):
            print(ps)
            print("{:.5f}".format(rts[n-1]))
    rts = clean_array(rts)
    return( rts)


def encontrar_divisor_quadratico(p_in):
    u=np.array([0.1,0.1])  # inicialização ps2_out= s^2+u(1)*s +u(2)
    du=np.zeros(2)
    for k in range (1,500):
       ps2_out= np.poly1d([1,u[0],u[1]])
       print(ps2_out)
       [pdiv, Rb, b] = fatorar_p2(p_in,ps2_out)
       pb_in=np.poly1d(np.flip(b))
       [pc   , Rc, c] = fatorar_p2(pb_in,ps2_out)
       du[0] =  ( b[0]*c[3]-c[2]*b[1] ) / (c[1]*c[3]-c[2]**2);
       du[1] =  ( c[1]*b[1]-b[0]*c[2] ) / (c[1]*c[3]-c[2]**2);
       print('erro=%.2e' % (np.linalg.norm(du)))
       if (np.linalg.norm(du)<1e-16):
           break
       u=u+du
    return (ps2_out,pdiv)


def bairstow_zeros(ps,prt):    
    n=len(ps)+1
    rts=np.zeros(n)*(1+1j)
    p_in=ps
    while (n>3): #repita enquanto ordem do polinômio > 2 
          [ps2_out,p_div]=encontrar_divisor_quadratico(p_in)
          (rts[n-1],rts[n-2],delta)=baskara(ps2_out);  
          if (prt):
              print(p_in)  
              print(ps2_out)
              print("{:.5f}".format(rts[n-2]))
              print("{:.5f}".format(rts[n-1]))
              print("\n")
          p_in=p_div
          n=n-2;   # diminua em 2 a ordem do polinômio
    if (prt):
        print(p_in)
    if(n==3):  # última equação ordem 2
        (rts[n-1],rts[n-2],delta)=baskara(p_in);
        if (prt):
            print("{:.5f}".format(rts[n-2]))
            print("{:.5f}".format(rts[n-1]))
    else:     # última equação ordem 1
         rts[n-1]=-p_in[0]/p_in[1]
         if (prt):
            print("{:.5f}".format(rts[n-1]))
    rts = clean_array(rts)
    return (rts)

