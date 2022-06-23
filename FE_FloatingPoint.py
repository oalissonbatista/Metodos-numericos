# https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
# https://matplotlib.org/3.5.0/tutorials/introductory/pyplot.html

import numpy as np
from tabulate import tabulate
import math
import matplotlib.pyplot as plt
from   FE_Inteiros import *

prt = True
n=8
m=24
e=np.array([0,1,0,0,0,1,0,0,1,1,0,1,0,0,1,0,1,0,0,1,0,1,0,1,1,1,0,0,0,0,1,0])
a=1684.679931640625

def dec2binario4(a,n,m): #real positivo 
  p= 2.0**np.arange(n-1,-m-1,-1)
  r=a/p
  r = r.astype(int)
  b = r%2
  return b

def binario2dec4(b,n,m): #real positivo 
   p= 2.0**np.arange(n-1,-m-1,-1)
   a= b@p  
   return a

def b2s4(b,n,m):
    return str(b[0:n])+'.'+str(b[n:n+m])

def float2dec(e,n,m): 
     bias = 2**(n-1)-1
     sinal= e[0] 
     E    = e[1:n+1]
     M    = np.append([1], e[n+1:n+m]) #  + 1 bit
     expoente10 = binario2dec1(E) - bias
     mantissa10 = binario2dec4(M,1,m-1)
     a = (-1)**sinal * (mantissa10) * 2.0**(expoente10)
     return a

def dec2float(a,n,m,prt):
    sinal = np.array([0]);
    if (a<0):
        sinal = np.array([1]);
    bias = 2**(n-1) -1;
    expoente10=-bias-1
    if (a!=0):
      expoente10 = np.floor(np.log2(np.abs(a)));
    if (expoente10<-bias):
        e= np.concatenate([sinal, np.zeros(n+m-1) ])
    elif (expoente10>bias+1): 
        e= np.concatenate([sinal, np.ones(n+m-1) ])
    else:
        E = dec2binario1(expoente10+bias,n);
        mantissa10= abs(a) * 2**(-expoente10);
        mantissa = dec2binario4(mantissa10,1,m-1);
        M = mantissa[1:m]  # descarta 1 bit

        e = np.concatenate([sinal, E, M])
    if (prt):
        printfloat(e,n,m)
    return e

def printfloat(e,n,m):
     print(str(e[0:1])+str(e[1:n+1])+'[1.]'+str(e[n+1:n+m]))

def printfloat_bitoculto(e,n,m,bit_oculto):
     if  (bit_oculto==0):
         prt_out=(str(e[0:1])+str(e[1:n+1])+'[0.]'+str(e[n+1:n+m]))
     elif (bit_oculto==1):
         prt_out=(str(e[0:1])+str(e[1:n+1])+'[1.]'+str(e[n+1:n+m]))
     else:
         prt_out=(str(e[0:1])+str(e[1:n+1])+'[(+1)0.]'+str(e[n+1:n+m]))
     return prt_out
    
def tabela_floating_point(n,m):
    Bias = 2**(n-1)-1
    E=np.zeros(n)
    tabela=np.zeros([2**(m-1)+1,2**n+1])
    tabela[0,0] = np.nan
    plus_one = np.array([0,1])
    E=np.zeros(n)
    for col in range(1,2**n+1):
        tabela[0,col] = binario2dec1(E)-Bias 
        E=SomaBinaria(E,plus_one,False)

    M=np.zeros(m-1)
    for lin in range (1,2**(m-1)+1):
        tabela[lin,0] = binario2dec4(np.concatenate([[1], M]),1,m-1)
        M=SomaBinaria(M,plus_one,False)

    E=np.zeros(n)
    for col in range(1,2**n+1):
        M=np.zeros(m-1)
        for lin in range (1,2**(m-1)+1):
            e=np.concatenate([[0], E, M])
            tabela[lin,col] = float2dec(e,n,m)
            M=SomaBinaria(M,plus_one,False)
        E=SomaBinaria(E,plus_one,False)
    pdtabulate=lambda tabela:tabulate(tabela,headers='keys')
    print(pdtabulate(tabela))  
    return tabela

def grafico_floating_point(n,m):
    a_in=np.zeros(2000)
    a_out=np.zeros(2000)
    max_x = float2dec(np.concatenate([[0],np.ones(n+m-1)]),n,m) 
    k=0
    for k in range (0,2000):
        x =  -1.5*max_x + 3.0*max_x * k/2000
        a_in[k] = x
        e=dec2float(a_in[k],n,m,False)
        a_out[k]=float2dec(e,n,m)
        k=k+1;
    plt.plot(a_in,a_out)
    plt.show()
    return

    # z = np.argwhere(y==1)
    # print(x[z[-1]+1]-x[z[0]])
def plot_epsilon_floating_point(n,m):
    x=np.linspace(0.90,1.3,10000)
    y=np.zeros(10000)
    for k in range (len(x)):
        e=dec2float(x[k],n,m,False);
        y[k]=float2dec(e,n,m)
    plt.plot(x,y)
    plt.show()
    return

def epsilon(n,m):
    dy = 1.0
    y = 2.0
    while (y!=1):
      dy=dy/2.0
      y=float2dec(dec2float(1.0+dy,n,m,False),n,m)
    eps=2.0*dy
    print(eps)
    return eps

def epsilon2(n,m):
    eps = 2**-(m-1)
    print(eps)
    return eps


def compatibiliza_entradas(x1,x2,prt):
   if(abs(x2)>abs(x1)): # Forçar abs(x1)>abs(x2) 
      temp=x1;
      x1=x2;
      x2=temp;
   trocar_sinais=0;
   if( x1<0 ):  #forçar x1>0
       trocar_sinais=1;
       x1=-x1;
       x2=-x2;
       if (prt):
         print("Trocar Sinais de x1 e x2\n" % (x1))
   return (x1,x2,trocar_sinais) 

def SomaFloat(x1_in,x2_in,n,m,prt): 
   (x1,x2,trocar_sinais)=compatibiliza_entradas(x1_in,x2_in,prt)
   e1=dec2float(x1,n,m,False) # e1 e e2 são os floats equivalente
   e2=dec2float(x2,n,m,False)
   exp1  = e1[1:n+1]    # expoentes
   exp2  = e2[1:n+1]
   mant1  = np.concatenate([np.array([0,1]), e1[n+1:n+m]]) # mantissas
   mant2  = np.concatenate([np.array([0,1]), e2[n+1:n+m]])
   desloc=binario2dec1(SomaBinaria(exp1,comple2(exp2),False))
   mant2d=np.zeros(m+1) #desloca mantissa direita para igualar expoentes
   mant2d[desloc:m+1]=mant2[0:m+1-desloc]
   e2_d = np.concatenate([[e2[0]],exp1, mant2d[2:m+1]])
   e2_d = e2_d.astype(int)
   if(prt): #imprime entradas e matissa deslocada
        if(abs(x2)<1e-5):
          print("(x1,x2)=(%e,%e)" % (x1,x2))
        else:
          print("(x1,x2)=(%f,%f)" % (x1,x2))
        str_out=printfloat_bitoculto(e1,n,m,1)
        print("e1     = %s" % (str_out))
        str_out=printfloat_bitoculto(e2,n,m,1)
        print("e2     = %s" % (str_out))
        if(desloc>=m):
          print("Underflow: deslocamento -> maior que %d" % (m)) 
        print("mantissa x2 deslocada -> de %d bits" % (desloc))
        str_out=printfloat_bitoculto(e2_d,n,m,mant2d[1])
        print("e2d    = %s" % (str_out))
   if (e2[0]==0):   # os dois números são positivos       
       mant_s = SomaBinaria(mant1,mant2d,False) #soma mantissas
       soma = np.concatenate([[0],exp1,mant_s[2:m+2]])
       if(mant_s[0]==1):   # vai 1 na soma
          soma_d = soma #deslocar mantissa direita e somar 1 ao expoente
          soma = np.concatenate([[0],SomaBinaria(exp1,[0,1],False),mant_s[1:m+1]])
       if(prt):
            str_out=printfloat_bitoculto(e2_d,n,m,mant2d[1])
            print("--------------------\ne2d    = %s" % str_out)
            str_out=printfloat_bitoculto(e1,n,m,1)
            print("e1     = %s" % str_out)
            print("-------------------- +")
            if (mant_s[0]==1):
              str_out=printfloat_bitoculto(soma_d,n,m,2)
              print("soma_d = %s" % str_out)
              print("Vai um! mantisa soma delocada de -> 1 bits")
   else: # o primeiro é positivo e o segundo número é negativo
       mant2d = mant2d.astype(int)
       mant_s = SomaBinaria(mant1,comple2(mant2d),False) #subtair mantissas
       soma_d = np.concatenate([[0], exp1, mant_s[2:m+2]])
     #  lista_uns = find(mant_s[1:m+1]==1)
       lista_uns = np.argwhere(mant_s[0:m+2]==1)
       if(lista_uns.size==0): # mantissa zerada
            soma = np.concatenate([[0], np.zeros(n+m-1)])
       else:   #deslocar mantissa para esquerda e subtrair expoente
           desloc=lista_uns[0] # encontrar primeiro '1'
           mant_sd=np.zeros([m+1])
           mant_sd=mant_sd.astype(int)
           desloc=desloc.astype(int)-1
           mant_sd[0:(m+2-desloc[0]-1)]=mant_s[desloc[0]:m+1]
           exp_s = SomaBinaria(exp1,dec2binario3(-desloc[0],n),False)
           soma = np.concatenate([[0], exp_s, mant_sd[2:m+1]])
       if (prt):
           mant2d_c2 = comple2(mant2d)
           temp = np.concatenate([[e1[1]],exp1,mant2d_c2[2:m+1]])
           str_out=printfloat_bitoculto(temp,n,m,mant2d_c2[1])
           print("--------------------\ne2d(c2)= %s" % (str_out))
           str_out=printfloat_bitoculto(e1,n,m,1)
           print("e1     = %s" % (str_out))
           str_out=printfloat_bitoculto(soma_d,n,m,mant_s[1])
           print("-------------------- +\nsoma_d = %s" % (str_out))
           if(lista_uns.size>0):
                print("mantissa soma deslocada de <- %d bits" % (desloc))
           else:
                print("Underflow: matissa nula após deslocamento <-")
   if (trocar_sinais):
     soma[0]=1
   soma = soma.astype(int)  
   x=float2dec(soma,n,m)
   if (prt):
       if (trocar_sinais):
         print("Trocar sinal da soma" % (x1))
       str_out=printfloat_bitoculto(soma,n,m,1)
       print("soma   = %s" % (str_out))
       print("Soma binária ( %d,%d)\t= %.12f" % (n,m,x)) 
       print("Soma exata   (11,53)\t= %.12f" % (x1_in+x2_in)) 
   return x
