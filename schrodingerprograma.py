
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 10:49:17 2023

@author: carlo
"""

from pylab import*
from numpy import*
import numpy as np
from numpy.linalg import*
import scipy.optimize as so
import cmath as cmath

import matplotlib.pyplot as plt
N=1000
n=0
nciclos=N/4
lamb=0.3
modulo=np.zeros(N)
conjugada=np.zeros(N,dtype=complex)
alpha=np.zeros(N,dtype=complex)
phi=np.zeros(N,dtype=complex)
amenos=np.ones(N,dtype=complex)
amas=np.ones(N,dtype=complex)
ao=np.zeros(N,dtype=complex)
V=np.zeros(N)
gamma=np.zeros(N,dtype=complex)
b=np.zeros(N,dtype=complex) 
beta=np.zeros(N,dtype=complex)
x=np.zeros(N,dtype=complex)
chi=np.zeros(N,dtype=complex)
norma=np.zeros(N)
r=np.zeros(N)
y=np.zeros(N)
ko=2*pi*nciclos/N
s=1/(4*ko**2) 


for i in range(N):
    phi[i]=np.exp(i*ko*1j)*np.exp((-8*(4*i-N)**2)/N**2)
phi[0]=0    
phi[N-1]=0


f=N-1

for i in range(N):
    if 2*N/5<=i<=3*N/5:
        V[i]= lamb*ko**2

    ao[i]=-2+(2/s)*1j-V[i]
i=N-1  
while  i>=0:
    gamma[i]=1/(ao[i]+amas[i]*alpha[i])
    alpha[i-1]=-amenos[i]*gamma[i]
    i=i-1
    
alpha[N-2]=0 #Rodrigo dice de poner N-2   
    
for i in range(len(phi)):
    phi[i]=np.exp(i*ko*1j)*np.exp((-8*(4*i-N)**2)/N**2)
    #b[i]=(4*phi[i]*1j)/s tiene que cambiarse para cada n pero si lo cambio me da nan
             
        
while n<=N: 
    for i in range(len(phi)):
        b[i]=(4*phi[i]*1j)/s #al ponerlo aquí phi se me va a inf
    
    i=N-1
    while i>=0:
         beta[i-1]=gamma[i]*(b[i]-amas[i]*beta[i])
         i=i-1
         
    beta[N-2]=0    
     
    
    for i in range(N-1):
        chi[i+1]=alpha[i]*chi[i]+beta[i]
        phi[i]=chi[i]-phi[i]   
    
    chi[N-1]=0
    chi[0]=0 
    phi[0]=0
    phi[N-1]=0
    for i in range(N):

        #modulo[i]=(phi.real[i]+phi.imag[i]*1j)*(phi.real[i]-phi.imag[i]*1j)
        modulo[i]=phi.real[i]**2+phi.imag[i]**2
    suma=0    
    for i in range(len(phi)):
        r[i]=phi[i].real
        y[i]=phi[i].imag
        norma[i]=(r[i]**2+y[i]**2)
        suma+=norma[i]
        
    norma1=suma**0.5      
    
    #g=np.arange(0, 1002, 1)
    with open("datosschrodingermodulo lambda=0.3 N=1000 n=250.txt", "a") as f:
        for i in range(len(phi)):
            f.write(f"{i}, {modulo[i]}, {V[i]}\n")   
            if i==N-1:
                f.write("\n") 
        f.write("\r\n")

          
            
            
    with open("normacte N=1000 n=250 lambda=0.3.txt", "a") as f:
        f.write(str(norma1))
        f.write(",")

    n=n+1

    #f=open("shrodinger.txt","w")
    #np.savetxt(f,n,delimiter=",")
    #np.savetxt(f,norma,delimiter=",")
    #np.savetxt(f,V, newline="\n")
    #data = np.column_stack((n, norma, V))
    #np.savetxt("shrodinger.txt", data, delimiter=",")
    # guardar: tiempo, la norma de phi, potencial
    #("%d, %d, %d", n, phi, V)
