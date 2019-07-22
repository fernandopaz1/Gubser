# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 09:09:29 2018

@author: Paz




Gusbser 

"""



import numpy as np
import scipy as sp
from scipy import special 
import math as mt
#from scipy import integrate as int
from scipy.integrate import quad, dblquad, nquad
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from mpmath import mp
from scipy.interpolate import InterpolatedUnivariateSpline 
from scipy.interpolate import InterpolatedUnivariateSpline as spline
#from goto import goto, label

plt.ion()

#from scipy.interpolate.InterpolatedUnivariateSpline import derivative as spinelder
#import sympy

#close("all")

t0=-10
tf=10
M=10000
N=2*M+1
h=(tf-t0)/N
tcambio=0
t1= np.linspace(t0,tcambio, N)
t2= np.linspace(tcambio+h, tf, N)
t = np.linspace(t0, tf, N)
t=np.concatenate((t1,t2),axis=0)
q2=1000
q=2*q2+1
indice1=q2

t_g=np.linspace(t0,tf,q)
              

h_g=(tf-t0)/q
zeta=np.zeros(N)
Temp=np.zeros(N)
zeta_lineal=np.zeros(N)


c=3.97887       # puede se c=39.7887 c=3.97887 c=1.1937 y c=0.397887
x0=1/np.sqrt(3)
T0=1.3468 
zeta0=0




q2=1000
q=2*q2+1
indice1=q2
tiemp=np.linspace(t0,tf,q)
              

hdnmr=(tf-t0)/q

integ_pi=np.zeros(N)
v_DNMR=np.zeros((2, q)) #v(1,:) es la energia y v(2,:) es el sh
T_DNMR=np.zeros(q)
T_DNMR_L=np.zeros(q)    



zeta[q2]=zeta0
Temp[q2]=T0  
zeta_lineal[q2]=zeta0



T_DNMR[0]=T0
#Esto es DNMR hago euler a partir del indice=indice0 (para adelante y para )
v_DNMR[0,q2]=(3/(np.pi**2))*(T0**4)
#esta ecuación esta sacada del paper Anisotropic fluid dynamics for a
#Gubser flow de Martinez Mcneils y Heinz ecuación 16
for i in range(0,q2-1):
        T_DNMR[q2+i]=(abs(((np.pi**2)/3)*v_DNMR[0,q2+i])**(1/4))
        v_DNMR[0,q2+i+1]=v_DNMR[0,q2+i]+ hdnmr*((4/3)*(v_DNMR[1,q2+i]-2)*v_DNMR[0,q2+i]*np.tanh(tiemp[q2+i]))
        v_DNMR[1,q2+i+1]=v_DNMR[1,q2+i]+ hdnmr*((-(T_DNMR[q2+i]/c)*(v_DNMR[1,q2+i]))+np.tanh(tiemp[q2+i])*((4/3)*(1/5+(5/14)*v_DNMR[1,q2+i]-((v_DNMR[1,q2+i])**2))))


for i in range(0,q2-1):
        T_DNMR[q2-i]=(abs(((np.pi**2)/3)*v_DNMR[0,q2-i])**(1/4))
        v_DNMR[0,q2-i-1]=v_DNMR[0,q2-i]- hdnmr*((4/3)*(v_DNMR[1,q2-i]-2)*v_DNMR[0,q2-i]*np.tanh(tiemp[q2-i]))
        v_DNMR[1,q2-i-1]=v_DNMR[1,q2-i]- hdnmr*((-(T_DNMR[q2-i]/c)*(v_DNMR[1,q2-i]))+np.tanh(tiemp[q2-i])*((4/3)*(1/5+(5/14)*v_DNMR[1,q2-i]-((v_DNMR[1,q2-i])**2)))) 


for i in range(0,q):
    T_DNMR_L[i]=(abs(((np.pi**2)/3)*v_DNMR[0,i])**(1/4))    
 






def H(x):
    if float(x)==1:
        a=1
    else:
        a=0.5*(x**2+(x**4)*((np.arctanh(np.sqrt(1-x**2))))/np.sqrt(1-x**2))
    return  a



def A(x):
    if float(x)==1:
        a=0
    else:
        a=(x*np.sqrt(x**2-1)*(1-2*x**2)+(1-4*x**2)*float(mp.acoth(x/np.sqrt(x**2-1))))/((2*x**3)*(x**2-1)**(3/2))
    return  a

H_aux=np.zeros(len(t))

for i in range(0,len(t)):
    H_aux[i]=H(np.cosh(0)/np.cosh(t[i]))



def D(t2,t1,spl_T):
    a=spl_T.integral(t1,t2)/c
    return np.exp(-a)


#T0=np.max(T_DNMR)

T=T0*np.ones(len(t_g))

T=T0*(np.cosh(0)/np.cosh(t_g))**(2/3)


spl_T=spline(t_g,T)


#D_aux=np.zeros(len(t_g))
#
#for i in range(0,len(t_g)):
#    D_aux[i]=D(t_g[i],-10,spl_T)

T1=np.zeros(len(t_g))
T2=np.zeros(len(t_g))

j=0
while j<1:
    for i in range(0,len(t_g)):
        T1[i]=D(t_g[i],0,spl_T)*H(np.cosh(0)/np.cosh(t_g[i]))*(T0**4)
        T2[i]=(1/c)*quad(lambda x: D(t_g[i],x,spl_T)*H(np.cosh(x)/np.cosh(t_g[i]))*(spl_T(x)**5), 0, t_g[i])[0]
        T[i]=(T1[i]+T2[i])**(1/4)
    spl_T=spline(t_g,T)
    j=j+1
        
#    
#plt.figure(1)
#plt.semilogy(t_g,T/T0,':',t_g,(np.cosh(0)/np.cosh(t_g))**(2/3),tiemp,T_DNMR_L/T_DNMR_L[q2])
#plt.show()


def z_fs(t,t_a,z_a):
    return (-1+(1+z_a)*(np.cosh(t_a)/np.cosh(t))**2)       

def R_200(x):
    if x<0:
        a=0.5*(1/(1+x)-(np.arctanh(-np.sqrt(abs(x)))/np.sqrt(abs(x))))
    elif x==0:
        a=1
    else:
        a=0.5*(1/(1+x)+(np.arctan(np.sqrt(x))/np.sqrt(x)))
    return a

def R_220(x):
    if x<0:
        a=0.5*(1/x)*(-1/(1+x)-(np.arctanh(-np.sqrt(abs(x)))/np.sqrt(abs(x))))
    elif x==0:
        a=1/3
    else:
        a=0.5*(1/x)*(-1/(1+x)+(np.arctan(np.sqrt(abs(x)))/np.sqrt(abs(x))))
    return a


#z_0=0
#x=t_g
#
#
#R_aux=np.zeros(len(t_g))
#for j in range(0,len(t_g)):
#    for i in range(0,len(t_g)):
#    #        R_aux[i]=quad(lambda x: D(t_g[i],x,spl_T)*((np.cosh(x)/np.cosh(t_g[i]))**4)*(spl_T(x))*E_rs(spl_T(x),z_fs(t_g[i],x,z_0)), 0, t_g[i])[0]
#        R_aux[i]=z_fs(t_g[i],x[j],z_0)
#    #        R_aux[i]=R_220(z_fs(t_g[i],x[j],z_0))
#    plt.plot(t_g,R_aux)


def E_rs(l,x):
    return (l**4)*R_200(x)

def pi_rs(l,x):
    return (l**4)*(R_220(x)-(1/3)*R_200(x))


z_0=0

lambda_0=T0

E1=np.zeros(len(t_g))
E2=np.zeros(len(t_g))
E=np.zeros(len(t_g))
pi1=np.zeros(len(t_g))
pi2=np.zeros(len(t_g))
pi_h=np.zeros(len(t_g))
pi_g=np.zeros(len(t_g))



#pi_aux=spl_T(t_g)*pi_rs(spl_T(t_g),z_fs(t_g,0,z_0))

#T0=0.002

#T=T0*np.ones(len(t_g))

#T=T0*(np.cosh(t0)/np.cosh(t_g))**(2/3)


spl_T=spline(t_g,T)
#
#j=0
#while j<100:
#    spl_T=spline(t_g,T)
#    for i in range(0,len(t_g)):
#            E1[i]=D(t_g[i],t0,spl_T)*((np.cosh(t0)/np.cosh(t_g[i]))**4)*E_rs(T0,z_fs(t_g[i],t0,z_0))
#            E2[i]=(1/c)*quad(lambda x: D(t_g[i],x,spl_T)*((np.cosh(x)/np.cosh(t_g[i]))**4)*(spl_T(x))*E_rs(spl_T(x),z_fs(t_g[i],x,z_0)), t0, t_g[i])[0]
#            T[i]=((E1[i]+E2[i]))**(1/4)
#    j=j+1
#

plt.semilogy(t_g,T/T0,tiemp,T_DNMR_L/T0,':')
#((np.pi**2)/3)

#plt.figure(2)
#plt.plot(t_g,E1,t_g,E2)
#plt.show()
#
#plt.figure(3)
#plt.semilogy(t_g,T,tiemp,T_DNMR_L,':')
#plt.show()
#               



for i in range(0,len(t_g)):
        pi1[i]=D(t_g[i],0,spl_T)*(A(np.cosh(t_g[i])/np.cosh(0)))*T0**4
        pi2[i]=(1/c)*quad(lambda x: D(t_g[i],x,spl_T)*(A(np.cosh(t_g[i])/np.cosh(x)))*(spl_T(x)**5), 0, t_g[i])[0]
        pi_h[i]=((pi1[i]+pi2[i]))/np.pi
        pi_g[i]=pi_h[i]/(T1[i]+T2[i])
j=j+1

#((np.pi**2)/3)

plt.figure(4)
plt.plot(t_g,(3/4)*pi_g,tiemp,v_DNMR[1,:],':')
plt.show()     


#plt.figure(5)
#plt.plot(t_g,pi1/(E1+E2),t_g,pi2/(E1+E2),':')
#plt.show()     
#
#plt.figure(6)
#plt.plot(t_g,pi1,t_g,pi2,':')
#plt.show()     

plt.figure(4)
plt.plot(t_g,((4/3)*(T1+T2)/T)*(np.cosh(t_g)**2))
plt.show() 
