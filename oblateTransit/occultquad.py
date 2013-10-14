from numpy import size,zeros,where,arccos,sqrt,pi,log
from cellip import *
#from ellip import *
import time
## tolerance for double precision equalities
## special case integrations
tol = 1e-14

def noplanet(nz):
    muo1 = zeros(nz) + 1. 
    mu0  = zeros(nz) + 1.
    return [muo1,mu0]

def oc_inegressi(z,p):

     x1=(p-z)**2.
     x2=(p+z)**2.
     tmp = (1.-p**2.+z**2.)/2./z
     tmp = where(tmp > 1.,1.,tmp)
     tmp = where(tmp < -1.,-1.,tmp)
     kap1 = arccos(tmp)
     tmp = (p**2.+z**2-1.)/2./p/z
     tmp = where(tmp > 1.,1.,tmp)
     tmp = where(tmp < -1.,-1.,tmp)
     kap0 = arccos(tmp)
     tmp = 4.*z**2-(1.+z**2-p**2)**2
     tmp = where(tmp < 0,0,tmp)
     lambdae = (p**2*kap0+kap1 - 0.5*sqrt(tmp))/pi
     # eta_1
     etad = 1./2./pi*(kap1+p**2*(p**2+2.*z**2)*kap0- \
           (1.+5.*p**2+z**2)/4.*sqrt((1.-x1)*(x2-1.)))
     return [lambdae,etad]
    
def oc_edge(z,p):
    if p < 0.5:
        # Case 5
        q=2.*p  # corrected typo in paper (2k -> 2p)
        Ek,Kk = ellke(q)
        # lambda_4
        lambdad = 1./3.+2./9./pi*(4.*(2.*p**2-1.)*Ek+\
                                       (1.-4.*p**2)*Kk)
        # eta_2
        etad = p**2/2.*(p**2+2.*z**2)        
        lambdae = p**2 # uniform disk
        return [lambdad,lambdae,etad]
    elif p > 0.5:
         # Case 7
        q=0.5/p # corrected typo in paper (1/2k -> 1/2p)
        Ek,Kk = ellke(q)
         # lambda_3
        lambdad = 1./3.+16.*p/9./pi*(2.*p**2-1.)*Ek-\
                       (32.*p**4-20.*p**2+3.)/9./pi/p*Kk
        lambdae =None
        etad = None
          # etad = eta_1 already
        return [lambdad,lambdae,etad]
    else:
         # Case 6
        lambdad = 1./3.-4./pi/9.
        etad = 3./32.
        lambdae = None
        return [lambdad,lambdae,etad]
def oc_iegress(z,p):

    # Case 2, Case 8 - ingress/egress (with limb darkening)
    x1=(p-z)**2.
    x2=(p+z)**2.
    x3=p**2.-z**2.
    q=sqrt((1.-x1)/(x2-x1))
    Ek,Kk = ellke(q)
    n=1./x1-1.

    # lambda_1:
    lambdad=2./9./pi/sqrt(x2-x1)* (((1.-x2)*(2.*x2+x1-3.)- \
            3.*x3*(x2-2.))*Kk+ (x2-x1)*(z**2+7.*p**2-4.)*Ek-\
            3.*x3/x1*ellpic_bulirsch(n,q))
    return lambdad

def oc_inside(z,p): 
    ## eta_2
    etad = p**2/2.*(p**2+2.*z**2)
  
    ## uniform disk
    lambdae = p**2
    return [lambdae,etad]

def oc_edgeedge(p):
    
    lambdad = 2./3./pi*arccos(1.-2.*p)-\
            4./9./pi*sqrt(p*(1.-p))*(3.+2.*p-8.*p**2)
    if p > 0.5:
        lambdad -= 2./3.
    return lambdad

def oc_between(z,p):

    x1=(p-z)**2.
    x2=(p+z)**2.
    x3=p**2.-z**2.
    q=sqrt((x2-x1)/(1.-x1))
    n=x2/x1-1.
    Ek,Kk = ellke(q)    
                
    ## Case 3, Case 9 - anywhere in between
    ## lambda_2
    lambdad = 2./9./pi/sqrt(1.-x1)*\
             ((1.-5.*z**2+p**2+x3**2)*Kk+\
             (1.-x1)*(z**2+7.*p**2-4.)*Ek-\
             3.*x3/x1*ellpic_bulirsch(n,q))
    return lambdad  

def fulloccult(z,u1,u2,p):
    #print notusedyet
    
    z = where(abs(p-z) < tol,p,z)
    z = where(abs((p-1)-z) < tol,p-1.,z)
 
    nz = size(z)
    lambdad = zeros(nz)
    etad = zeros(nz)
    lambdae = zeros(nz)
    omega=1.-u1/3.-u2/6.

 #   x1=(p-z)**2.
 #   x2=(p+z)**2.
 #   x3=p**2.-z**2.

    indexoccult = z <= p-1. #,complement=notused2)
    occult = z[indexoccult]
    iengressi = z[-indexoccult]
    if indexoccult.any():
        etad[indexoccult] = 0.5 
        lambdae[indexoccult] = 1.0
    if indexoccult.all():
        muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.* \
                                             (p > z))+u2*etad)/omega
        mu0=1.-lambdae
        return [muo1,mu0]
    else:
        # Case 2,7,8 - ingress,egress, uniform disk only
        lambdae[-indexoccult],etad[-indexoccult]= oc_inegressi(iengressi,p)
        # Case 5, 6, 7 - the edge of planet lies at origin of star
        ocltor = z == p #, complement=notused3)
        if (ocltor.any()): 
            lambdad_oc,lambdae_oc,etad_oc = oc_edge(z[ocltor],p)
            lambdad[ocltor] = lambdad_oc
            if (not etad_oc == None):
                etad[ocltor]= etad_oc
            if (not lambdae_oc == None):
                lambdad[ocltor]= lambdad_oc
        if(ocltor.all()):
            muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*\
                       (lambdad+2./3.*(p > z))+u2*etad)/omega
            mu0=1.-lambdae
            return [muo1,mu0]
        else:
        # Case 2, Case 8 - ingress/egress (with limb darkening)
            #print lambdae,etad
            indexinegress = (-ocltor)*(-indexoccult)
            lambdad[indexinegress] = oc_iegress(z[indexinegress],p)
            #print len(z[indexinegress])
            muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.*\
                    (p > z))+u2*etad)/omega
            mu0=1.-lambdae
            #print muo1,mu0
            return [muo1,mu0]

def fulldisk(z,u1,u2,p):
  
    nz = size(z)
    lambdad = zeros(nz)
    etad = zeros(nz)
    lambdae = zeros(nz)
    omega=1.-u1/3.-u2/6.
     
    z = where(abs(p-z) < tol,p,z)
    z = where(abs((1-p)-z) < tol,1.-p,z)
    z = where(z < tol,0.,z)
    #print time.time()  
    indexoccult = z < 1-p  #,complement=notused2)
    iengressi = z[-indexoccult]
    if size(iengressi) != 0:
        # Case 2, 7, 8 - ingress/egress (uniform disk only)
        lambdae[-indexoccult],etad[-indexoccult]= oc_inegressi(iengressi,p)
   
    # Case 5, 6, 7 - the edge of planet lies at origin of star
    ocltor = z == p #, complement=notused3)

    if (ocltor.any()):  
        lambdad_oc,lambdae_oc,etad_oc = oc_edge(z[ocltor],p)
        lambdad[ocltor] = lambdad_oc
        if (not etad_oc == None):
            etad[ocltor]= etad_oc
        if (not lambdae_oc == None):
            lambdad[ocltor]= lambdad_oc
    if(ocltor.all()):
        muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*\
                   (lambdad+2./3.*(p > z))+u2*etad)/omega
        mu0=1.-lambdae
        return [muo1,mu0]
    else:
        #print time.time()  
        ## Case 4 - edge of planet hits edge of star
        indexedge = z == 1.-p #, complement=notused6)
        if indexedge.any():
         ## lambda_5
            lambdad[indexedge] = oc_edgeedge(p)
        if indexedge.all():
             muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*\
                       (lambdad+2./3.*(p > z))+u2*etad)/omega
             mu0=1.-lambdae
             return [muo1,mu0]
        else: 
    # Case 2, Case 8 - ingress/egress (with limb darkening)
            indexinegress = (-ocltor)*(-indexoccult)*(-indexedge)
            lambdad[indexinegress] = oc_iegress(z[indexinegress],p)
            indexinside = indexoccult + (p<0.5)*(z>=p) 
            if (not indexinside.any()):
                muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.*\
                                                          (p > z))+u2*etad)/omega
                mu0=1.-lambdae
                return [muo1,mu0]
            else:
                #print "reach case 3,4,9,10"
                # Case 3, 4, 9, 10 - planet completely inside star
                #print time.time()  
                if indexoccult.any():
                    lambdae[indexoccult],etad[indexoccult]=oc_inside(z[indexoccult],p)
                
                ## Case 10 - origin of planet hits origin of star
                origin = z == 0#, complement=notused7)
                if origin.any():
                    ## lambda_6
                    lambdad[origin] = -2./3.*(1.-p**2)**1.5
                if origin.all():
                    muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*\
                            (lambdad+2./3.*(p > z))+u2*etad)/omega
                    mu0=1.-lambdae
                    return [muo1,mu0]
                else:
                    #anywhereinbetween
                    indexbetween = indexoccult *(-indexedge)*(-origin)*(-ocltor)
                    lambdad[indexbetween] = oc_between(z[indexbetween],p) 
                #    print lambdad
                    p0=p #!CH hack here    
                    if p0 > 0:
                        muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.*\
                                                                (p > z))+u2*etad)/omega
                        mu0=1.-lambdae
                    else:
                        muo1 =1.+((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.*\
                                                                (p > z))+u2*etad)/omega
                        mu0=1.+lambdae
                    return [muo1,mu0]

# Computes Hasting's polynomial approximation for the complete
# elliptic integral of the first (ek) and second (kk) kind
#   Python translation of IDL code.
#   This routine computes the lightcurve for occultation of a
#   quadratically limb-darkened source without microlensing.  Please
#   cite Mandel & Agol (2002) and Eastman et al (2012) if you make use
#   of this routine in your research.  Please report errors or bugs to
#   jeastman@lcogt.net
def occultquad(z,u1,u2,p0):
    #print time.time()
    nz = size(z)


    p = abs(p0)
   ## trivial case of no planet
    if p <= 0.:
        return noplanet(nz)
    ## Case 1 - the star is unocculted:
    ## only consider points with z lt 1+p
    #notusedyet = where( z < (1. + p) )
    index = z<(1.+p)
    mu1outtran,mu0outtran = occultouttran(z[-index],u1,u2,p)
    #print '1',time.time()
    if p >= 1.:
        mu1intran,mu0intran = fulloccult(z[index],u1,u2,p)
    else:
        mu1intran,mu0intran = fulldisk(z[index],u1,u2,p)

    #print '2',time.time()
    #mu1intran,mu0intran = occultintran(z[index],u1,u2,p)
    mu1 = zeros(nz)
    mu0 = zeros(nz)
    mu1[index] = mu1intran
    mu1[-index] = mu1outtran
    mu0[index] = mu0intran
    mu0[-index] = mu0outtran
    #print time.time()
    return [mu1,mu0]

def occultouttran(z,u1,u2,p):
        nz = size(z)
        lambdad = zeros(nz)
        etad = zeros(nz)
        lambdae = zeros(nz)
        omega=1.-u1/3.-u2/6.
        muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.*(p > z))+ \
                  u2*etad)/omega
        mu0=1.-lambdae
        return [muo1,mu0]


def occultintran(z,u1,u2,p):
    notusedyet = where( z < (1. + p) )
    notusedyet = notusedyet[0]
    nz = size(z)

    # Case 11 - the  source is completely occulted:

