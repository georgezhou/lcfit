import os
import sys
import functions
from fitting_functions import *
from numpy import *
import matplotlib.pyplot as plt
from scipy import optimize
import random

def fake_transit(planet_f_i=0.0, planet_alpha_i=45.0, bsigma=1.*10**(-4), gsigma=1.*10**(-4), ntransits=2):

    ### Global constants
    au = 1.496*10**11
    msun = 1.988435*10**30
    rsun = 6.955*10**8
    mjup = 1.8988*10**27 
    rjup = 6.9173*10**7
    day = 60.*60.*24.
    gconst = 6.67*10**(-11)

    ### Model parameters
    ld1_coeff = eval(functions.read_config_file("LC_LD1"))
    ld2_coeff = eval(functions.read_config_file("LC_LD2"))
    period_i = float(functions.read_config_file("PERIOD"))
    t0_i = float(functions.read_config_file("T0"))
    rsum_i = float(functions.read_config_file("RSUM"))
    rratio_i = float(functions.read_config_file("RRATIO"))
    i_0_i = float(functions.read_config_file("I0"))
    ecosw_i = float(functions.read_config_file("ECOSW"))
    esinw_i = float(functions.read_config_file("ESINW"))
    edepth_i = float(functions.read_config_file("EDEPTH"))
    beta_i = float(functions.read_config_file("BETA"))
    fratio_i = float(functions.read_config_file("FRATIO"))
    theta_i = float(functions.read_config_file("THETA"))
    phi_i = float(functions.read_config_file("PHI"))
    Protot_i = float(functions.read_config_file("ROTP"))
    #planet_f_i = float(functions.read_config_file("PLANET_F"))
    #planet_alpha_i = float(functions.read_config_file("PLANET_ALPHA"))

    input_params = [ld1_coeff,ld2_coeff,period_i,t0_i,rsum_i,rratio_i,i_0_i,ecosw_i,esinw_i,edepth_i,beta_i,fratio_i,theta_i,phi_i,Protot_i,planet_f_i,planet_alpha_i]

    mstar = float(functions.read_config_file("MSTAR"))*msun
    rstar = 1.0*rsun

    cadence = "short"

    tdur = (period_i / pi) * arcsin((rsum_i**2+cos(i_0_i*pi/180.)**2)**0.5)
    
    hjd_list = []
    for i in range(ntransits):
        hjd_i = arange(t0_i+period_i*i-1.*tdur,t0_i+period_i*i+1.*tdur,1.0/(24.*60.))
        hjd_list.append(hjd_i)

    #################################

    from math import sqrt
    from scipy.stats import norm

    def brownian(x0, input_array, delta):

        n = len(input_array)
        dt = (input_array[-1]-input_array[0])/float(n)

        x0 = asarray(x0)

        # For each element of x0, generate a sample of n numbers from a
        # normal distribution.
        r = norm.rvs(size=x0.shape + (n,), scale=delta*sqrt(dt))

        out = empty(r.shape)

        # This computes the Brownian motion by forming the cumulative sum of
        # the random samples. 
        cumsum(r, axis=-1, out=out)

        # Add the initial condition.
        out += expand_dims(x0, axis=-1)

        return out

    def gaussian(input_array,mu,sigma):
        o = zeros(len(input_array))
        for i in range(len(o)):
            o[i] = random.gauss(mu,sigma)
        return o

    def make_model(hjd_i,input_params):
        ld1_coeff,ld2_coeff,period_i,t0_i,rsum_i,rratio_i,i_0_i,ecosw_i,esinw_i,edepth_i,beta_i,fratio_i,theta_i,phi_i,Protot_i,planet_f_i,planet_alpha_i = input_params

        rmeanrstar = rratio_i
        rratio_i = sqrt(1-planet_f_i)*rmeanrstar

        sma = rsum_i/(1+rmeanrstar)
        sma = 1/sma

        rsum_mean = rsum_i
        rsum_i = (1/sma) * (1+rratio_i)
        #print rsum_i

        ### Input for transit model 
        # ! V(1) = surface brightness ratio    V(15) = third light
        # ! V(2) = sum of fractional radii     V(16) = phase correction
        # ! V(3) = ratio of stellar radii      V(17) = light scaling factor
        # ! V(4) = linear LD for star A        V(18) = integration ring size (deg)
        # ! V(5) = linear LD for star B        V(19) = orbital period (days)
        # ! V(6) = orbital inclination         V(20) = ephemeris timebase (days)
        # ! V(7) = e cos(omega) OR ecentricity V(21) = nonlinear LD for star A
        # ! V(8) = e sin(omega) OR omega       V(22) = nonlinear LD for star B
        # ! V(9) = gravity darkening 1         
        # ! V(10)= gravity darkening 2         
        # ! V(11) = primary reflected light    
        # ! V(12) = secondary reflected light  
        # ! V(13) = stellar mass ratio        
        # ! V(14) = tidal lead/lag angle (deg)

        model_input = array([edepth_i,rsum_i,rratio_i,ld1_coeff[0],0,i_0_i,ecosw_i,esinw_i,0,0,0,0,0,0,0,0,0,1,period_i,t0_i,ld2_coeff[0],0])

        #print len(model_input),"model_input_length"
        model = transitmodel(model_input,hjd_i,1.,0.,cadence)
        #plt.plot(hjd_i,model)

        ### Apply obliquity
        rplanet = rratio_i*rstar
        a = (rplanet+rstar)/rsum_i
        b = a * cos(i_0_i*pi/180.)

        phase = (hjd_i-t0_i)/period_i
        phase = phase - floor(phase)

        if min(phase) < 0.05 or max(phase)>0.95:

            if cadence == "short":
                obliq_model = obliquity(hjd_i,t0_i,beta_i,fratio_i,theta_i,phi_i,Protot_i,b,a,period_i,rstar,mstar)
            else:
                obliq_model = []
                for datapoint in hjd_i:
                    datapoint = arange(datapoint-0.0104,datapoint+0.0104,0.00208)
                    #datapoint = arange(datapoint-0.0104,datapoint+0.0104,0.0005)

                    obliq_model.append(mean(obliquity(datapoint,t0_i,beta_i,fratio_i,theta_i,phi_i,Protot_i,b,a,period_i,rstar,mstar)))

                obliq_model = array(obliq_model)

            model = (model-max(model))*(1-fratio_i)*obliq_model+1

        #plt.plot(hjd_i,model)
        #plt.show()

        #sys.exit()


        ### Apply planet oblation

        model = model - oblateness_func.oblateness_func(hjd_i,t0_i,period_i,rmeanrstar,planet_f_i,planet_alpha_i,sma,i_0_i,ld1_coeff[0],ld2_coeff[0])

        return model


    lc = []
    for hjd_i in hjd_list:

        model = make_model(hjd_i,input_params)

        bnoise = brownian(0., hjd_i, bsigma)
        gnoise = gaussian(hjd_i,0.,gsigma)

        #model = model + bnoise + gnoise

        model += gnoise

        data = transpose([hjd_i,model,ones(len(hjd_i))])

        # plt.scatter(hjd_i,model,s=0.1)
        # plt.show()

        lc.append(data)

    return lc

if __name__ == "__main__":
    
    fake_transit()
