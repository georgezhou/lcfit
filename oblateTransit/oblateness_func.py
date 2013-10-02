import numpy as np
import scipy as sp
from math import sqrt,cos
import oblateness as Obl
import matplotlib
from matplotlib import pyplot as plt
def oblateness_func(hjd,t0,period,rmean,f,alpha,sma,inc,u1,u2):
    # rmean = 0.1
    # f = 0.1
    # alpha = 1/180.*np.pi
    # sma = 8.924 ### a/rstar
    # period = 2.218573 
    # inc = 90/180.*np.pi
    # u1 = 0.076
    # u2 = 0.034
    # Npoint = 500
    # percent = 0.1

    phase = (hjd-t0)/period
    phase = phase - np.floor(phase)
    ind = phase > 0.5
    phase[ind] -= 1.0

    alpha = alpha * np.pi / 180.
    inc = inc * np.pi / 180.

    Npoint = len(phase)
    dflux = np.zeros(Npoint)
    req = rmean/sqrt(1-f)
    rpol = sqrt(1-f)*rmean
    #initial
    obl = Obl.Oblateness(req,rpol,alpha,sma,inc,u1,u2)
    b0 = sma*cos(inc)
    if(b0>(1+req)):
        print 'no transit happen'
        return
    #dphi = 2*percent/(Npoint-1)
    totalFlux = np.pi*(1.0-u1/3-u2/6)
    #phi=-1*percent+dphi*np.arange(Npoint)
    #call
    obl.relativeFlux(phase,dflux)
    #plt.plot(phase,dflux/totalFlux)
    #plt.show()
    return dflux/totalFlux
