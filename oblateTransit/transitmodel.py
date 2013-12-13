import numpy as np
from occultquad import occultquad
from math import sqrt,cos
import matplotlib.pyplot as plt

def transitmodel(model_params,hjd,dummy1,dummy2,cadence):

    ### model_params t0 period rratio rsum inc u1 u2

    #t0,period,rpol,rsum,inc,u1,u2 = model_params
    edepth_i,rsum,rpol,u1,dummy,inc,ecosw_i,esinw_i,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,period,t0,u2,dummy=model_params

    sma = rsum/(1+rpol)
    sma = 1/sma
    inc = inc * np.pi/180.
    
    if cadence == "short":
        phi = (hjd - t0)/period

        z=sma*np.sqrt(np.sin(phi*2*np.pi)*np.sin(phi*2*np.pi)+np.cos(phi*2*np.pi)*np.cos(phi*2*np.pi)*cos(inc)*cos(inc));
        circularflux = occultquad(z,u1,u2,rpol)[0]
        index = np.cos(phi*2*np.pi)<0
        circularflux[index]=1.0
        flux = circularflux


    if cadence == "long":
        flux = []
        for i in hjd:
            delt = 0.00208
            phi = np.arange(i-0.0104,i+0.0104+delt,delt)
            phi_orig = phi
            phi = (phi-t0)/period

            z=sma*np.sqrt(np.sin(phi*2*np.pi)*np.sin(phi*2*np.pi)+np.cos(phi*2*np.pi)*np.cos(phi*2*np.pi)*cos(inc)*cos(inc));
            circularflux = occultquad(z,u1,u2,rpol)[0]
            index = np.cos(phi*2*np.pi)<0
            circularflux[index]=1.0

            flux.append(np.mean(circularflux))
        flux = np.array(flux)

    return flux

if __name__ == "__main__":
    model_params = [0.0,1.0,0.1,0.05,90.,0.25,0.3]
    hjd = np.arange(-0.05,0.05,0.002)

    flux = transitmodel(hjd,model_params,"long")
    plt.plot(hjd,flux)
    plt.show()
