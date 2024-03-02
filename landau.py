# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 18:12:06 2023
LANDAU LEVELS
@author: Alvaro Cauqui Diaz
"""
import numpy as np
import matplotlib.pyplot as plt

hbar=1.054571817e-34
pi=np.pi
m=9.109e-31
e=1.6e-19
B=hbar/e
w=(e*B)/m
xi0=(m*w)/hbar

lb=np.sqrt(hbar/e*B)
def xi(pos):
    return xi0*pos
def phi0(xi_):
    return (xi0/pi)**(0.25)*np.exp(((-xi_**2)/2*xi0))
def phi1(xi_):
    return (4/xi0*pi)**(0.25)*xi_*np.exp(((-xi_**2)/2*xi0))

X=np.arange(-20,20,0.00001)
xii=xi(X)




wavefunctions=np.zeros((100,len(X)))
wavefunctions[0]=phi0(xii)
wavefunctions[1]=phi1(xii)
for i in range(2,100,1):
    wavefunctions[i]=np.sqrt(2/i)*(xii*wavefunctions[i-1]-np.sqrt((i-1)/2)*wavefunctions[i-2])
#%%
plt.figure()
plt.plot(xii,wavefunctions[17])
plt.title('n=17')
plt.ylabel('$\phi_{17}$')
plt.xlabel('position')
plt.figure()
plt.plot(xii,wavefunctions[40])
plt.title('n=40')
plt.ylabel('$\phi_{40}$')
plt.xlabel('position')
plt.figure()
plt.plot(xii,wavefunctions[99])
plt.title('n=99')
plt.ylabel('$\phi_{99}$')
plt.xlabel('position')
plt.figure()
plt.plot(xii,wavefunctions[17],label='n=17')
plt.plot(xii,wavefunctions[40],label='n=40')
plt.plot(xii,wavefunctions[99],label='n=99')
plt.legend()

#%%
plt.figure()
plt.plot(xii,wavefunctions[17]*wavefunctions[17])
plt.title('n=17')
plt.ylabel('$|\phi_{17}|^2$')
plt.xlabel('position')
plt.figure()
plt.plot(xii,wavefunctions[40]*wavefunctions[40])
plt.title('n=40')
plt.ylabel('$|\phi_{40}|^2$')
plt.xlabel('position')
plt.figure()
plt.plot(xii,wavefunctions[99]*wavefunctions[99])
plt.title('n=99')
plt.ylabel('$|\phi_{99}|^2$')
plt.xlabel('position')
#%%
def classical(A,x):
    return (1/pi)*(1/(np.sqrt(A**2-x**2)))
X0=np.arange(-15,15,0.0000075)
classic=[classical(21,d) for d in X0]
A=1/sum(wavefunctions[99])
plt.figure()
plt.plot(xii,classic)
#plt.plot(xii,A*wavefunctions[99]*wavefunctions[99])

