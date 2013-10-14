#!/usr/bin/env python
from elliptic import Elliptic 
import numpy as np
elp = Elliptic()
def ellke(k):
    if(not type(k)==float):
        ee = np.zeros(len(k))*1.0
        kk = np.zeros(len(k))*1.0
        elp.Calkk(k,kk)
        elp.Calee(k,ee)
    else:
        k = np.array([k])
        ee = np.array([0.0])
        kk = np.array([0.0])
        elp.Calee(k,ee)
        elp.Calkk(k,kk)
        ee = ee[0]
        kk = kk[0]
    return [ee,kk]

def ellpic_bulirsch(n,k):
    if(not type(n)==float):
        nk = np.zeros(len(n))
        elp.Ellpic_bulirsch(n,k,nk)
        return nk
    else:
        n = np.array([n])
        k = np.array([k])
        nk = np.array([0.0])
        elp.Ellpic_bulirsch(n,k,nk)
        return nk[0]

