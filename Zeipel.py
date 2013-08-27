#!/usr/bin/env python
import numpy as np
import scipy as sp
import math
import matplotlib 
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib import cm

def derterminant(x,y,f,phi,Req):
    """x,y is the cooridnates refer to star center
    f is the oblatness
    phi is the stellar obliquity, assume phi is already in the right 
    units here
    Req is the equator radius of the star
    """
    da = 4.*y**2.*(1.-(1.-f**2))**2*np.sin(phi)**2.*np.cos(phi)**2.
    db = np.cos(phi)**2*(1-f)**2.+np.sin(phi)**2.
    dc = (y**2*np.sin(phi)**2-Req**2+x**2)*(1-f)**2+y**2*np.cos(phi)**2.
    #print 'ds',da,db,dc
    d = da-4*db*dc
    #print d
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

def cal_geff(x,y,f,phi,Req,ggraveq,groteq):
    """
    the effective gravity that counts for temperature change is a 
    combination of both the gravity and the rotation.
    """
    d=derterminant(x,y,f,phi,Req)
    z=cal_zcoord(x,y,f,phi,d)
    #print d,z
    x0,y0,z0=cal_rotated_coord(x,y,z,phi)
    #print 'new coords',x0,y0,z0
    absR = np.sqrt(x0**2.+y0**2.+z0**2.)
    absRper = np.sqrt(x0**2.+z0**2.)
    #print 'absR=',absR,'absRper=',absRper
    gi = -ggraveq*(Req/absR)**2*x0/absR + groteq/(Req/absRper)*x0/absRper
    gj = -ggraveq*(Req/absR)**2.*y0/absR
    gz = -ggraveq*(Req/absR)**2.*z0/absR + groteq/(Req/absRper)*z0/absRper
    g = np.sqrt(gi**2+gj**2+gz**2)
    return g


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

def Zeipel():
    #set up the test case

    M = 1.8; Req0 = 2.029 ; #M is in solar mass and R is is solar radius 
    f = 0.1947
    Tpole = 8450; beta = 0.1947; 
    Protot = 8.2 #Protot is in hour
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
    ellipse = Ellipse([0,0],2*Req,2*Rpole,fc='none',ec='k')
    for phi in phis:
        phi /=(180./math.pi)
        g = cal_Zeipel(X,Y,gratio,f,phi)
        T = Tpole*(g/gpole)**beta
        index = 170+count
        count+=1
        ax = fig.add_subplot(index)
        exterior = np.sqrt(((X/Req)**2) + ((Y/Rpole)**2)) > 1.0
        T_mask = np.ma.masked_array(T,mask=exterior)
        CS = ax.contourf(X*Req0,Y*Req0,T,cmap=cm.copper)
        cbar = plt.colorbar(CS)
        ax.set_xlim([-3,3])
        ax.set_ylim([-3,3])
        ax.add_patch(ellipse)
    plt.show()
     
    return

if __name__=='__main__':
    Zeipel()
