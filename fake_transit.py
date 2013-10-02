import os
import sys
import functions
from fitting_functions import *
from numpy import *
import matplotlib.pyplot as plt
from scipy import optimize
import random

### Global constants
au = 1.496*10**11
msun = 1.988435*10**30
rsun = 6.955*10**8
mjup = 1.8988*10**27 
rjup = 6.9173*10**7
day = 60.*60.*24.
gconst = 6.67*10**(-11)

### Model parameters
ld1_coeff = [0.25]
ld2_coeff = [0.30]
period_i = 100.0
t0_i = 0.0
rsum_i = 0.02
rratio_i = 0.1
i_0_i = 89.5
ecosw_i = 0.0
esinw_i = 0.0
edepth_i = 0.0
beta_i = 0.0
fratio_i = 0.0
theta_i = 0.0001
phi_i = 0.0
Protot_i = 100000.
planet_f_i = 0.3
planet_alpha_i = 45.

mstar = 1.0*msun
rstar = 1.0*rsun


cadence = "short"
hjd_i = arange(-0.6,0.6,1.0/(24.*60.))



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

def make_model():
    global ld1_coeff,ld2_coeff,period_i,t0_i,rsum_i,rratio_i,i_0_i,ecosw_i,esinw_i,edepth_i,beta_i,fratio_i,theta_i,phi_i,Protot_i,planet_f_i,planet_alpha_i

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

    
model = make_model()

bnoise = brownian(0., hjd_i, 1*10.**(-4))
gnoise = gaussian(hjd_i,0.,1*10.**(-4))

model = model + bnoise + gnoise

data = transpose([hjd_i,model,ones(len(hjd_i))])

o = open("data/transit_"+str(planet_f_i)+"_"+str(planet_alpha_i),"w")
functions.write_table(data,o)
o.close()


plt.scatter(hjd_i,model,s=0.1)
plt.show()
