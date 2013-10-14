#!/usr/bin/env python
import math
import numpy as np
import scipy as sp
import os
from dataio import *
from HATlc import lightcurve as lc
from HATlc import tran
import random
import oblateness as Obl
from occultquad import occultquad
#from Eastman_occultquad import occultquad
import matplotlib
from matplotlib import pyplot as plt
import cmd_parse as cmdp
import cfg_parse as cfgp
def binlc(filltime,time,mag):
    cadence = 30./60./24.
    nt = len(filltime)
    fillmag = np.zeros(nt)
    for i in xrange(nt):
        samples=mag[np.abs(time-filltime[i])<cadence/2.]
        if(len(samples)>3):
            fillmag[i] = np.mean(samples)
    return fillmag

def gentran(time,period,epoch,q):
    ftime=sp.zeros(len(time))
    ftime=(time-epoch-0.5*period)/period-((time-epoch-0.5*period)/period).astype(int)
    ind=ftime<0
    ftime[ind]+=1
    intran=(ftime > (0.5-q/2.0))*(ftime < (0.5+q/2.0))
    return intran


def trangen(time,mag,transit,lcflag=False):
    '''functions to generate fake lightcurves'''
    rmean = math.sqrt(transit.dip)
    f = transit.f
    inc = transit.inc*math.pi/180.
    alpha = transit.alpha*np.pi/180.
    u1 = transit.u1
    u2 = transit.u2
    sma = 1./transit.q/np.pi
    #print sma,transit.P,transit.u1,transit.u2,transit.epoch,transit.q,transit.qg
    req = rmean/math.sqrt(1-f)
    rpol = math.sqrt(1-f)*rmean
    ntran=int((max(time)-transit.epoch)/transit.P)-int((min(time)-transit.epoch)/transit.P)+1
    
       
    tranwidth=transit.q*transit.P/2.
    obl = Obl.Oblateness(req,rpol,alpha,sma,inc,u1,u2)
    totalFlux = np.pi*(1.0-u1/3-u2/6)
    intran = gentran(time,transit.P,transit.epoch,transit.q)
    medianmag=np.median(mag[-intran])
    if transit.stdmag==-1:
        stdmag=np.std(mag[-intran])
    else:
        stdmag=transit.stdmag
    #print rpol,medianmag
    phase = (time-transit.epoch)/transit.P-((time-transit.epoch)/transit.P).astype(int)

    ind = phase>0.5
    phase[ind]-=1.0
    fkmag=medianmag+np.random.randn(len(time))*stdmag
    if(lcflag):
        #initial short cadence data
        circularfluxmeanlc=medianmag+np.random.randn(len(time))*stdmag
        cadence = 1./60./24.
        nt = (int)((max(time)-min(time))/cadence)
        fktime = min(time)+np.arange(nt)*cadence
        fkmagsc=medianmag+np.random.randn(len(fktime))*stdmag
        
        #initial short cadence pahse
        fkintran = gentran(fktime,transit.P,transit.epoch,transit.q) 
        fkphase = abs((fktime-transit.epoch)/transit.P)-(abs((fktime-transit.epoch)/transit.P)).astype(int)
        ind = fkphase>0.5
        fkphase[ind]-=1.0
        z=sma*np.sqrt(np.sin(fkphase*2*np.pi)*np.sin(fkphase*2*np.pi)+np.cos(fkphase*2*np.pi)*np.cos(fkphase*2*np.pi)*np.cos(inc)*np.cos(inc))
        
        #compute corrections
        fkdflux = np.zeros(len(fktime[fkintran]))
        obl.relativeFlux(fkphase[fkintran],fkdflux)
        #fkdflux/=totalFlux
        circularfluxmean = occultquad(z,u1,u2,rmean)[0]
        circularfluxmean = 1-(1-circularfluxmean)/(totalFlux/np.pi)
        circularflux = occultquad(z,u1,u2,rpol)[0]
        circularflux = 1-(1-circularflux)/(totalFlux/np.pi)
        index = np.cos(fkphase*2*np.pi)<0
        circularflux[index]=1.0
        circularfluxmean[index]=1.0
        
        fkmagsc[fkintran] = medianmag
        fkmagsc[fkintran]+=1-(circularflux[fkintran]-fkdflux)
        circularfluxmean=1-circularfluxmean+medianmag
        fkmag[intran] = binlc(time[intran],fktime[fkintran],fkmagsc[fkintran])
        residual = (fkmagsc-circularfluxmean)/1.e-6
        circularfluxmeanlc[intran] = binlc(time[intran],fktime[fkintran],circularfluxmean[fkintran]) + np.random.randn(len(time[intran]))*stdmag
        circularfluxmean = circularfluxmeanlc 
        plt.plot(fkphase[fkintran],residual[fkintran],'r') 
        plt.xlim([-0.01,0.01])

    else:
        dflux = np.zeros(len(time[intran]))
        obl.relativeFlux(phase[intran],dflux)
        dflux/=totalFlux
        z=sma*np.sqrt(np.sin(phase*2*np.pi)*np.sin(phase*2*np.pi)+np.cos(phase*2*np.pi)*np.cos(phase*2*np.pi)*np.cos(inc)*np.cos(inc))
    
        circularflux = occultquad(z,u1,u2,rpol)[0]
    
        index = np.cos(phase*2*np.pi)<0
        circularflux[index]=1.0
    
        circularfluxmean,check = occultquad(z,u1,u2,rmean)
        circularfluxmean[index]=1.0
        
        plt.plot(phase[intran],(circularfluxmean[intran]-circularflux[intran]+dflux)/1.e-6,'r')
        plt.xlim([-0.01,0.01])
        fkmag[intran] = medianmag+np.random.randn(len(fkmag[intran]))*stdmag
        fkmag[intran]+=1-(circularflux[intran]-dflux)
        circularfluxmean=1-circularfluxmean+medianmag+np.random.randn(len(fkmag))*stdmag

    plt.plot(phase[intran],(-circularfluxmean[intran]+fkmag[intran])/1.e-6,'o',mec='b',mfc='None',ms=1.5,mew=1)
    plt.show()
    
    return [fkmag,circularfluxmean]

	
def main():
    lcfile=lc()
    transit=tran()
    options = cmdp.ltf_parse()
    infileflag = options.infileflag
    outfileflag = options.outfileflag
    inpathflag = options.inpathflag
    noplot=options.noplot
    cfgfile = options.cfg
    uflag = options.uflag
    lcflag = int(cfgp.File_parse(cfgfile,'lcflag'))
    infile = cfgp.File_parse(cfgfile,'infile')
    inpath = cfgp.File_parse(cfgfile,'inpath')
    outfile = cfgp.File_parse(cfgfile,'outfile')
    coljd= int(cfgp.File_parse(cfgfile,'coljd'))
    colmag= int(cfgp.File_parse(cfgfile,'colmag'))
    if not inpath=='':
        os.chdir(inpath)
    transit.readpara(cfgfile)
    time=[];readcolumn(time,coljd,infile);time=np.array(time)
    mag=[];readcolumn(mag,colmag,infile);mag=np.array(mag)
    fkmag,circularflux=trangen(time,mag,transit,lcflag=lcflag)
    if outfile=='':
        outfile=os.path.splitext(infile)[0]+'.fktrap'
    fout=open(outfile,mode='w')
    fout.write('#Generated by lcgen with the following paramters:\n')
    fout.write('#Period = %f, rpstar=%f, Tdur= %f, inc=%f, alpha=%f, f=%f\n' % (transit.P,math.sqrt(transit.dip),transit.q*transit.P*24.,transit.inc,transit.alpha,transit.f))
    for i in range(len(time)):
        fout.write('%13.6f %13.6f %13.6f %13.6f\n' % (time[i],mag[i],fkmag[i],circularflux[i]))
    fout.close()
    return


if __name__=='__main__':
	main()
