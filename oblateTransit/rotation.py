#!/usr/bin/env python
#first attempt to convert f to Prot
from math import pi,sqrt
import numpy as np
import scipy as sp
import matplotlib 
from matplotlib import pyplot as plt
from const import *
def SeagerHui(f,Req,Mp):
    #assuming Req and Mp already in cgs.
    Prot = 2*pi*np.sqrt(Req**3./(2*f*G*Mp))
    return Prot

def CarterWinn(f,Req,Mp,J2):
    Prot = 2*pi*np.sqrt(Req**3./(G*Mp*(2*f-3*J2)))
    return Prot

def DarwinRadau(f,C=0.23):
    return (-0.3+2.5*C-15./8.*C**2.)*f

def Prot_cal(f,Req,Mp,J2flag=False):
    if (J2flag):
        J2 = DarwinRadau(f)
        return CarterWinn(f,Req,Mp,J2)
    else:
        return SeagerHui(f,Req,Mp)

def main():
    solarf = np.array([0.00012,0.00009,0.0035,0.0052,0.06487,0.09796,0.02293,0.01708])
    J2 = np.array([0.00006,0.000004,0.001083,0.00196,0.014736,0.016298,0.003343,0.003411])
    Prot_real = np.array([58.6462,243.0187,0.99726968,1.025986,0.41007,0.426,0.71833,0.67125])*day
    Req = np.array([0.383,0.95,1.0,0.532,10.97,9.14,3.98,3.86])*rearth
    Mp = np.array([330.2e24,4868.5e24,5973.6e24,641.85e24,1.8986e30,568460.0e24,86832.0e24,102430.0e24])
    Prot1 = CarterWinn(solarf,Req,Mp,J2)/day
    Prot2 = Prot_cal(solarf,Req,Mp,J2flag=True)/day
    Prot3 = Prot_cal(solarf,Req,Mp)/day
    xgrid = np.arange(len(solarf))+1
    plt.semilogy(xgrid,Prot_real/day,label='observed')
    plt.semilogy(xgrid,Prot1,label='measured J2')
    plt.semilogy(xgrid,Prot2,label='fitted J2')
    plt.semilogy(xgrid,Prot3,label='no J2')
    plt.xlabel('Sequence in solar')
    plt.ylabel('Protation(day)')
    plt.legend()
    plt.show()

    return

if __name__=='__main__':
    main()
