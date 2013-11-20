#!/usr/bin/env python
import numpy as np
import scipy as sp
#from scipy.optimize import broyden1
#from scipy.optimize import newton_krylov
#from scipy.optimize import root
from scipy.optimize import fsolve
import math
import matplotlib 
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib import cm

def cal_tau(theta,w,r):
    tau = 1./3.*w**2.*r**3.*np.cos(theta)**3.+np.cos(theta)+np.log(np.tan(theta/2.))
    #tau = np.cos(theta)+np.log(np.tan(theta/2.))
    return tau

def solveR(theta,w,Rabs):
    def F(x,theta,w):
        return 1./(w**2.*x)+0.5*x**2*np.sin(theta)**2.-1./w**2.-0.5
    def DF(x,theta,w):
        return -1./(w**2.*x**2.)+x*np.sin(theta)**2.
    r = np.zeros(theta.shape)
    for i in range(len(r)):
        r[i]= fsolve(F,Rabs[i],args=(theta[i],w),fprime=DF)
    
    return r

def solveangle(tau,theta):
    def F(x,tau0):
        return (1-x**2.)/(1+x**2.)+np.log(x)-tau0
    
    def DF(x,tau0):
        return -4*x/(1+x**2.)**2.+1./x
    
    guess = np.tan(theta/2.)
    
    x = np.zeros(tau.shape)
    iflag = np.zeros(tau.shape)
    for i in range(len(x)):
        x[i]= fsolve(F,guess[i],args=tau[i],fprime=DF)
   
    return x


def cal_Tratio(theta,r,f,w):
    Tpole = (1-f)**(-0.5)*(np.exp(2./3.*w**2.*(1-f)**3.))**(0.25)
    tau = cal_tau(theta,w,r)
    
    Thalf = solveangle(np.ravel(tau),np.ravel(theta)).reshape(tau.shape) 
    
    tanTheta = (2*Thalf)/(1-Thalf**2.)
    indexa = abs(np.tan(theta))<1.e-7
    T = np.zeros(theta.shape)
    T[indexa]= 1/r[indexa]**0.5*(np.exp(2./3.*w**2.*r[indexa]**3.))**0.25
    indexb = abs(theta-np.pi/2.)<1.e-7
    T[indexb] = (1/r[indexb]**4.+w**4.*r[indexb]**2-2*w**2./r[indexb])**(1./8.)*(1-w**2.*r[indexb]**3.)**(-1./6.)
    Ta = np.zeros(theta.shape)
    Tb = np.zeros(theta.shape)
    #print np.where(r==0)
    Ta[-indexa*-indexb] = 1/r[-indexa*-indexb]**4.
    Ta[-indexa*-indexb]+=w**4.* r[-indexa*-indexb]**2.*np.sin(theta[-indexa*-indexb])**2.
    Tb[-indexa*-indexb] = 2*w**2.*np.sin(theta[-indexa*-indexb])**2./r[-indexa*-indexb]
    Ta[-indexa*-indexb]-= Tb[-indexa*-indexb]
    T[-indexa*-indexb]=Ta[-indexa*-indexb]**(1./8.)
    T[-indexa*-indexb]*=np.sqrt(tanTheta[-indexa*-indexb]/np.tan(theta[-indexa*-indexb]))
    
    return T/Tpole

def derterminant(x,y,f,phi,Req):
    """x,y is the cooridnates refer to star center
    f is the oblatness
    phi is the stellar obliquity, assume phi is already in the right 
    units here
    Req is the equator radius of the star
    The power of f in Barnes 2009 is wrong
    """
    da = 4.*y**2.*(1.-(1.-f)**2.)**2*np.sin(phi)**2.*np.cos(phi)**2.
    db = np.cos(phi)**2*(1-f)**2.+np.sin(phi)**2.
    dc = (y**2*np.sin(phi)**2-Req**2+x**2)*(1-f)**2+y**2*np.cos(phi)**2.
    d = da-4*db*dc
    return d

def cal_zcoord(x,y,f,phi,d):
    """
    compute the vertical coordinate z.
    """
    za = -2*y*(1-(1-f)**2)*np.sin(phi)*np.cos(phi)+np.sqrt(d)
    zb = 2*((1-f)**2*np.cos(phi)**2.+np.sin(phi)**2.)
    z = za/zb

    return z

def cal_rotated_coord(x,y,z,phi):
    """
    x0,y0,z0 is the rotated coordinate from x,y,z in the y-z plan so that 
    y0 (the rotation axises) is inside the plane perpenticular to line of sight. 
    """
    x0 = x
    y0 = y*np.cos(phi)+z*np.sin(phi)
    z0 = -y*np.sin(phi)+z*np.cos(phi)
    return [x0,y0,z0]

def cal_Lara(x,y,f,phi,groteq,Req):
    """
    the effective gravity that counts for temperature change is a 
    combination of both the gravity and the rotation.
    """
    #print x,y,f,phi,Req
    d=derterminant(x,y,f,phi,Req)
    z=cal_zcoord(x,y,f,phi,d)
    exterior = np.isnan(z)
    x_mask = np.ma.masked_array(x,mask=exterior)
    y_mask = np.ma.masked_array(y,mask=exterior)
    d_mask = np.ma.masked_array(d,mask=exterior)
    z_mask=np.ma.masked_array(z,mask=exterior)
     
    x0,y0,z0=cal_rotated_coord(x,y,z,phi)
    
    absR = np.sqrt(x0**2.+y0**2.+z0**2.)
    absRper = np.sqrt(x0**2.+z0**2.)
    #print 'absR=',absR,'absRper=',absRper
    theta = np.arccos(y0/absR) 
    
    w = np.sqrt(groteq) 
    #tempR = solveR(np.ravel(theta),w,np.ravel(absR)).reshape(theta.shape)

    #print theta,w
    ratio = cal_Tratio(theta,absR,f,w)
    #ratio= cal_Tratio(theta,tempR,f,w)
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #CS = ax.contourf(x,y,ratio,cmap=cm.copper)
    #ax.set_title('phi=%f' % phi)
    #ax.set_xlim([-1.5,1.5])
    #ax.set_ylim([-1,1])
    #plt.show()
    return ratio

def cal_Zeipel(X,Y,gratio,f,phi):
    """
    The function to connect with other routines.
    X,Y, the coordinates that normalized to Req; 
    f, the oblateness of the star; 
    phi, the obliquity of the star spin to the plain perpenticular to 
    line of sight;
    gratio, groteq/ggraveq, the ratio of the effective accelation due 
    to rotation versus gravity accelaration at equator. 
    """
    g = cal_geff(X,Y,f,phi,1.0,1.0,gratio)
    return g

def Lara():
    #set up the test case

    M = 1.8; Req0 = 2.029 ; #M is in solar mass and R is is solar radius 
    f = 0.1947
    Tpole = 8450; beta = 0.1947; 
    Protot = 8.64 #Protot is in hour
    G_MRsun_hour = 4*math.pi**2*(1.5e3/7)**3./(365.*24.)**2
    Omega = 2*math.pi/Protot
    
    Rpole0 = Req0*(1-f) 
    ggraveq0 = G_MRsun_hour*M/Req0**2.
    groteq0 = Omega**2.*Req0
    ggraveq = 1.0
    groteq = ggraveq/ggraveq0*groteq0
    gratio = groteq/ggraveq
    
    Req = 1.0
    Rpole =1-f 
    x = 0; y = Rpole
    gpole=ggraveq*Req**2/Rpole**2.
    x = np.linspace(-Req,Req,100)
    y = np.linspace(-Rpole,Rpole,100)
    X, Y = np.meshgrid(x, y)
    phis = [0.0,30.0,60.0,90.0,120.0,150.0,180.0]
    
    fig = plt.figure(figsize=[35,3.5])
    count = 1
    ellipse = Ellipse([0,0],2*Req0,2*Rpole0,fc='none',ec='k')
    for phi in phis:
        phi /=(180./math.pi)
        ratio = cal_Lara(X,Y,f,phi,gratio,Req)
        T = Tpole*ratio
        index = 170+count
        count+=1
        ax = fig.add_subplot(index)
        #ax = fig.add_subplot(111)
        exterior = np.sqrt(((X/Req)**2) + ((Y/Rpole)**2)) > 1.0
        T_mask = np.ma.masked_array(T,mask=exterior)
        CS = ax.contourf(X*Req0,Y*Req0,T,cmap=cm.copper)
        cbar = plt.colorbar(CS)
        ax.set_xlim([-3,3])
        ax.set_ylim([-3,3])
        ax.add_patch(ellipse)
    #plt.tight_layout()
    plt.show()
     
    return

if __name__=='__main__':
    Lara()
