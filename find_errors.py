import os
import sys
import matplotlib.pyplot as plt
from numpy import *
import functions
from scipy import interpolate
import string

Gconst = 6.674*10**(-8)
Msun = 1.988*10**33
Rsun = 6.955*10**10


### Get free parameters:


### Read from config file
### Load initial parameters and lightcurve

temp_param_names = []
temp_param_vals = []
temp_param_range = []

#lc = functions.read_config_file("INPUT_LC")
#lc = loadtxt(lc)

lc_ld1 = eval(functions.read_config_file("LC_LD1"))
lc_ld1_err = eval(functions.read_config_file("LC_LD1_ERR"))
lc_ld2 = eval(functions.read_config_file("LC_LD2"))
lc_ld2_err = eval(functions.read_config_file("LC_LD2_ERR"))
for i in range(len(lc_ld1)):
    temp_param_names.append("lc_ld1")
    temp_param_vals.append(lc_ld1[i])
    temp_param_range.append(lc_ld1_err[i])
    temp_param_names.append("lc_ld2")
    temp_param_vals.append(lc_ld2[i])
    temp_param_range.append(lc_ld2_err[i])

temp_param_names.append("period")
temp_param_vals.append(float(functions.read_config_file("PERIOD")))
temp_param_range.append(float(functions.read_config_file("PERIOD_ERR")))

temp_param_names.append("t0")
temp_param_vals.append(float(functions.read_config_file("T0"))-floor(float(functions.read_config_file("T0"))))
temp_param_range.append(float(functions.read_config_file("T0_ERR")))

temp_param_names.append("rsum")
temp_param_vals.append(float(functions.read_config_file("RSUM")))
temp_param_range.append(float(functions.read_config_file("RSUM_ERR")))

temp_param_names.append("rratio")
temp_param_vals.append(float(functions.read_config_file("RRATIO")))
temp_param_range.append(float(functions.read_config_file("RRATIO_ERR")))

temp_param_names.append("i_0")
temp_param_vals.append(float(functions.read_config_file("I0")))
temp_param_range.append(float(functions.read_config_file("I0_ERR")))

temp_param_names.append("ecosw")
temp_param_vals.append(float(functions.read_config_file("ECOSW")))
temp_param_range.append(float(functions.read_config_file("ECOSW_ERR")))

temp_param_names.append("esinw")
temp_param_vals.append(float(functions.read_config_file("ESINW")))
temp_param_range.append(float(functions.read_config_file("ESINW_ERR")))

temp_param_names.append("edepth")
temp_param_vals.append(float(functions.read_config_file("EDEPTH")))
temp_param_range.append(float(functions.read_config_file("EDEPTH_ERR")))

temp_param_names.append("beta")
temp_param_vals.append(float(functions.read_config_file("BETA")))
temp_param_range.append(float(functions.read_config_file("BETA_ERR")))

temp_param_names.append("fratio")
temp_param_vals.append(float(functions.read_config_file("FRATIO")))
temp_param_range.append(float(functions.read_config_file("FRATIO_ERR")))

temp_param_names.append("theta")
temp_param_vals.append(float(functions.read_config_file("THETA")))
temp_param_range.append(float(functions.read_config_file("THETA_ERR")))

temp_param_names.append("phi")
temp_param_vals.append(float(functions.read_config_file("PHI")))
temp_param_range.append(float(functions.read_config_file("PHI_ERR")))

temp_param_names.append("Protot")
temp_param_vals.append(float(functions.read_config_file("ROTP")))
temp_param_range.append(float(functions.read_config_file("ROTP_ERR")))

temp_param_names.append("planet_f")
temp_param_vals.append(float(functions.read_config_file("PLANET_F")))
temp_param_range.append(float(functions.read_config_file("PLANET_F_ERR")))

temp_param_names.append("planet_alpha")
temp_param_vals.append(float(functions.read_config_file("PLANET_ALPHA")))
temp_param_range.append(float(functions.read_config_file("PLANET_ALPHA_ERR")))


temp_param_names.append("tdiff")
temp_param_vals.append(float(functions.read_config_file("TDIFF")))
temp_param_range.append(float(functions.read_config_file("TDIFF_ERR")))


### Distribute as free or fixed param, including associated brackets
free_param_names = []
free_param_vals = []
free_param_range = []
free_param_func = []

fixed_param_names = []
fixed_param_vals = []

for i in range(len(temp_param_names)):
    if temp_param_range[i] == 0:
        fixed_param_names.append(temp_param_names[i])
        fixed_param_vals.append(temp_param_vals[i])
    else:
        free_param_names.append(temp_param_names[i])
        free_param_vals.append(temp_param_vals[i])
        free_param_range.append(temp_param_range[i])
        free_param_func.append("b")


axes_tested = free_param_names

tested = loadtxt("mcmc_tested_params")
prob = loadtxt("mcmc_chisq_log")
prob = transpose(prob)

minlen = min([len(tested),len(prob)])
tested = tested[:minlen]
prob = prob[:minlen]

prob = exp((min(prob)-prob)/2)

best_param = []
maxprob = max(prob)
for i in range(len(prob)):
    if prob[i] == maxprob:
        best_param = tested[i]
        break
print "best parameters",best_param

tested = transpose(tested)


def bin(nbins,xlist,zlist):
    xmin = min(xlist)
    xmax = max(xlist)

    x_range = xmax-xmin
    xstep = x_range/nbins
    
    xaxis = arange(xmin,xmax+xstep/2,xstep)

    binarray = []
    for i in range(len(xaxis)):
        binarray.append([])

    for i in range(len(zlist)):
        xbin = int((xlist[i]-xmin)/xstep)
        binarray[xbin].append(zlist[i])

    zarray = []
    for i in range(len(binarray)):
        if len(binarray[i]) > 10:
            sortarray = sort(binarray[i])
            sortarray = sortarray[-50:]
            value = median(sortarray)
            if functions.isnan(value):
                zarray.append(value)
            else:
                zarray.append(0)
            
            # zarray.append(max(binarray[i]))
        else:
            zarray.append(0)

    zarray = array(zarray)
    zarray = zarray / zarray.max()
    xaxis = xaxis + xstep/2
    return xaxis,zarray

def find_edges(x,y):
    xmin = min(x)
    xmax = max(x)
    for i in range(len(x)):
        if y[i] > 0.606:
            xmin = x[i]
            break

    for i in range(len(x)):
        if y[i] > 0.606:
            xmax = x[i]
    return xmin,xmax

def fitsmooth(x,y,box):
    #print box
    ytemp = []
    for i in range(len(x)):
        imin = i-box
        imax = i+box+1
        if imin <0:
            imin = 0
        if imax >len(x):
            imax = len(x)
        ybox = y[imin:imax]
        ytemp.append(mean(ybox))

    f = interpolate.splrep(x,ytemp)
    x = arange(min(x),max(x),(max(x)-min(x))/1000)
    f = interpolate.splev(x,f)
    f = f /max(f)

    return x,f

for i in range(len(axes_tested)):
    try:
        plt.clf()
        #plt.scatter(tested[i],prob)
        x,y = bin(500,tested[i],prob)
        plt.plot(x,y,"r-")    

        if axes_tested[i] in ["ecosw","esinw","rratio","rsum","t0","k0","period","q0","i_0"]:
            x,fit = fitsmooth(x,y,5)
        else:
            x,fit = fitsmooth(x,y,10)

        xmin,xmax = find_edges(x,fit)
        plt.axvline(x=xmin)
        plt.axvline(x=xmax)

        xmid = x[0]
        for j in range(len(x)):
            if fit[j] == max(fit):
                xmid = x[j]
                break

        #xmid = best_param[i]

        plt.axvline(x=xmid)

        print axes_tested[i],best_param[i],xmid,xmin-xmid,xmax-xmid
        #best_param.append(xmid)
        best_param[i] = xmid

        plt.plot(x,fit,"k-")    
        plt.xlabel(axes_tested[i])
        plt.show()
    except ValueError:
        pass

best_param = [axes_tested,best_param]
f = open("best_param_mcmc","w")
functions.write_table(best_param,f)
f.close()

