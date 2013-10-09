#!/usr/peyton/common/software/python/2.7.2/bin/python

import os,sys,functions
from numpy import *
import fake_transit,fit_lc


def fit_fake_transit(planet_f=0.0, planet_alpha=45.0):

    ### Define noise and ntransit
    brownian_noise=1.5*10**(-4)
    gaussian_noise=1.5*10**(-4)
    no_transits=2

    ### Make transit
    lc = fake_transit.fake_transit(planet_f_i=planet_f,planet_alpha_i=planet_alpha,bsigma=brownian_noise,gsigma=gaussian_noise,ntransits=no_transits)

    ### Fit the transit and return xmin,xmid,xmax for planet_f
    xmin,xmid,xmax = fit_lc.wrap(lc)

    o = open("output","a")
    o.write(str(planet_f)+" "+str(xmin)+" "+str(xmid)+" "+str(xmax)+"\n")
    o.close()

    print planet_f,xmin,xmid,xmax

def start_process(config_file):
    config_file = loadtxt(config_file)
    print config_file
    fit_fake_transit(planet_f=config_file[0],planet_alpha=config_file[1])

if __name__ == "__main__":
    #fit_fake_transit()
    start_process(sys.argv[1])
