# https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
# https://matplotlib.org/3.5.0/tutorials/introductory/pyplot.html

import numpy as np
import matplotlib.pyplot as plt
f = lambda x: 4.9 + (5.1*x-22.3)/(3.1*x-4.2)
f = lambda x: 4*np.exp(x/3)-20*np.exp(-x/5)*np.sin(x)
f = lambda x: x**2-8*x+6
f = lambda x: x**5-6*x**4+10*x*3-10*x**2+5*x-2
#print(np.array([x,f(x)]).transpose())


def Bissecao_abs(f,a,b,tol,prt): #erro absoluto #f=função a,b sao intervalos / tol
    if (f(a)*f(b)>0):
       print("não há raiz no intervalo [%f,%f]]" % (a,b))
       return []
    if (prt): 
        print('i \ta\t\t\tx1\t\t\tb\t\t\terro') 
    for k in range(1,200):   
        erro=abs(b-a)
        x1=(a+b)/2     #Bisseção
        if (prt): 
            print("%d\t%.10f(%2d)\t%.10f(%2d)\t%.10f(%2d)\t%.1e"
            % (k,a,np.sign(f(a)),x1,np.sign(f(x1)),b,np.sign(f(b)),erro) )
        if  ( (erro<tol) | (f(x1)==0) ):
            break
        if  (f(x1)*f(a)<0):
            b=x1
        else:
            a=x1
    return x1

def Bissecao_n(f,a,b,tol,prt): #número de iterações
    if (f(a)*f(b)>0):
       print("não há raiz no intervalo [%f,%f]]" % (a,b))
       return []
    n=np.ceil(np.log2((b-a)/tol))+1
    if (prt): 
        print("Previsão de %d iterações" % (n)) 
        print ('i \ta\t\t\tx1\t\t\tb\t\t\terro')
    for k in range(1,int(n+1)): 
        x1=(a+b)/2    # Bisseção
        if (prt): 
            print("%d\t%.10f(%2d)\t%.10f(%2d)\t%.10f(%2d)\t%.1e"
            % (k,a,np.sign(f(a)),x1,np.sign(f(x1)),b,np.sign(f(b)),abs(b-a)) ) 
        if  ( f(x1)==0 ):
            break
        if  (f(x1)*f(a)<0):
            b=x1
        else:
            a=x1         
    return x1

def Bissecao(f,a,b,tol,prt): #erro relativo
    if (f(a)*f(b)>0):
       print("não há raiz no intervalo [%f,%f]]" % (a,b))
       return []
    if (prt): 
        print ('i \ta\t\t\tx1\t\t\tb\t\t\terro') 
    x1=a;
    erro=  1;
    for k in range(1,200):  
        x0=x1
        x1=(a+b)/2  # Bisseção
        if(x1!=0):
            erro =abs((x1-x0)/x1)
        if (prt): 
            print ( '%i\t%.10f(%2d)\t%.10f(%2d)\t%.10f(%2d)\t%.1e'
            % (k,a,np.sign(f(a)),x1,np.sign(f(x1)),b,np.sign(f(b)),erro) )
        if  ( (erro<tol) | (f(x1)==0) ):
            break;
        if  (f(x1)*f(a)<0):
            b=x1
        else:
            a=x1
    return x1


def FalsaPosicao(f,a,b,tol,prt):
    if (f(a)*f(b)>0):
       print("não há raiz no intervalo [%f,%f]]" % (a,b))
       return []
    if (prt): 
        print ('i \ta\t\t\tx1\t\t\tb\t\t\terro')
    x1=a;
    erro=  1;
    for k in range(1,500): 
        x0=x1
        x1=(a*f(b)-b*f(a))/(f(b)-f(a)) #Falsa Posição (Bisseção x1=(a+b)/2)
        if(x1!=0):
            erro =abs((x1-x0)/x1)
        if (prt): 
            print ( '%i\t%.10f(%2d)\t%.10f(%2d)\t%.10f(%2d)\t%.1e'
            % (k,a,np.sign(f(a)),x1,np.sign(f(x1)),b,np.sign(f(b)),erro) )
        if  ( (erro<tol) | (f(x1)==0) ):
            break
        if  (f(x1)*f(a)<0):
            b=x1
        else:
            a=x1
    return x1



def FalsaPosicaoModificado(f,a,b,tol,prt):
    if (f(a)*f(b)>0):
       print("não há raiz no intervalo [%f,%f]]" % (a,b))
       return []   
    if (prt): 
        print ('i \ta\t\t\tx1\t\t\tb\t\t\terro') 
    x1=a; 
    fa = f(a);
    fb = f(b); 
    na=0; 
    nb=0; 
    for k in range(1,500): 
        x0=x1;
        x1=(a*fb-b*fa)/(fb-fa) # falsa posição modificado
        fx1 = f(x1)
        if(x1!=0):
            erro =abs((x1-x0)/x1)
        if (prt): 
            print ( '%i\t%.10f(%2d)\t%.10f(%2d)\t%.10f(%2d)\t%.1e'
            % (k,a,np.sign(f(a)),x1,np.sign(f(x1)),b,np.sign(f(b)),erro) )
        if  ( (erro<tol) | (f(x1)==0) ):
            break
        if (fx1*fa < 0):
            b=x1
            fb = fx1;
            nb=0;
            na=na+1
            if (na>=2):
                fa = fa/2
        else:
            a=x1;
            fa = fx1;
            na=0;
            nb=nb+1
            if (nb>=2):
                fb = fb/2
    return x1

f = lambda x: x**2-8*x+6
a=0
b=1
tol = 0.0001
prt = True
