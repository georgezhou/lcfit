import os,sys,functions
from numpy import *
import fake_transit,fit_lc


def fit_fake_transit(planet_f=0.0, planet_alpha=45.0):

    ### Define noise and ntransit
    brownian_noise=1.*10**(-4)
    gaussian_noise=1.*10**(-4)
    no_transits=2

    ### Make transit
    lc = fake_transit.fake_transit(planet_f_i=planet_f,planet_alpha_i=planet_alpha,bsigma=brownian_noise,gsigma=gaussian_noise,ntransits=no_transits)

    ### Fit the transit and return xmin,xmid,xmax for planet_f
    xmin,xmid,xmax = fit_lc.wrap(lc)

    return planet_f,xmin,xmid,xmax

if __name__ == "__main__":
    fit_fake_transit()
