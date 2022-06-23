# https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
# https://matplotlib.org/3.5.0/tutorials/introductory/pyplot.html

import numpy as np
import math
import matplotlib.pyplot as plt

def dec2binario1(a,n): 
    p= 2**np.arange(n-1,-1,-1)
    r=a/p
    r = r.astype(int)
    b = r%2 
    if((a>(2**n)-1)|(a<0)):
        b=[]
    return b

def binario2dec1(b): 
    p= 2**np.arange(b.size-1,-1,-1);
    a= b@p  
    return a

def  dec2binario2(a,n): 
  b=dec2binario1(abs(a),n-1)
  if (a>=0):
      b=np.concatenate( ([0], b ) ) # acrescentar o bit de sinal
  else:
      b=np.concatenate( ([1], b ) )
  if(  a<(-2^(n-1)-1) |  a>(2^(n-1)-1)  ):
      b=[]
  return b

def binario2dec2(b):
   sinal=b[0]
   a= (-1)**sinal*binario2dec1(b[1::]);  
   return a


def Compatibiliza_Tamanhos(b1,b2):
    n1=len(b1)
    n2=len(b2)
    b11=b1
    b12=b2
    if (n1>n2): #length(b1)<>length(b2)
        b2=np.concatenate( (np.ones(n1-n2)*b2[0], b2) )
        b2 = b2.astype(int)
    elif (n2>n1):
        b1=np.concatenate( (np.ones(n2-n1)*b1(0), b1) )
        b1 = b1.astype(int)
    return (b1,b2)



def SomaBinaria(b1,b2,prt):
    (b1,b2)=Compatibiliza_Tamanhos(b1,b2)
    n=len(b1)
    v1=0
    b3=np.zeros(n)
    for k in range(len(b1)-1,-1,-1):
       b3[k]=b1[k]+b2[k]+v1
       v1=0
       if (b3[k]==2):
           b3[k]=0
           v1=1
       elif (b3[k]==3):
           b3[k]=1
           v1=1
       b3 = b3.astype(int)
    if (prt):
        print("  %s\n+ %s\n= %s \n" % (str(b1),str(b2),str(b3)))
    return b3

def comple2(b):
    b=~b+2   # inverter bits
    b=SomaBinaria(b,[0, 1],False)  # somar 1    
    return b
 
def dec2binario3(a,n): 
    b=dec2binario1(abs(a),n)
    if(a<0):
        b=comple2(b)
    if(  a<(-2^(n-1)) |  a>(2^(n-1)-1)  ):
        b=[] 
    return b 

def binario2dec3(b):
    sinal = b[0]
    if(sinal):
        b=comple2(b)
    a= (-1)**sinal * binario2dec1(b) ;
    return a


prt=True
b1=np.array([1, 1, 0, 0 ,1, 0, 1, 0, 1,1])
b2=np.array([1, 1, 1 ,0, 0, 1, 0])
