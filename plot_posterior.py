import os
import sys
import functions
from numpy import *
import matplotlib.pyplot as plt
from scipy import interpolate


def bin1d(nbins,xlist,zlist,createaxis=True,xaxis=0):

    if createaxis:

        xlist_temp = []
        zlist_temp = []
        for i in range(len(xlist)):
            if zlist[i] > 0.001:
                xlist_temp.append(xlist[i])
                zlist_temp.append(zlist[i])
        xlist,zlist = xlist_temp,zlist_temp
    
        xmin = min(xlist)
        xmax = max(xlist)

        x_range = xmax-xmin
        xstep = x_range/nbins
    
        xaxis = arange(xmin,xmax+xstep/2,xstep)

    else:
        xmin = min(xaxis)
        xmax = max(xaxis)
        x_range = xmax-xmin
        xstep = x_range/nbins   

    binarray = []
    for i in range(len(xaxis)):
        binarray.append([])

    for i in range(len(zlist)):
        xbin = int((xlist[i]-xmin)/xstep)
        if xbin < len(binarray)-1 and xbin >= 0:
            binarray[xbin].append(zlist[i])

    zarray = []
    for i in range(len(binarray)):
        if len(binarray[i]) > 10:
            sortarray = sort(binarray[i])
            sortarray = sortarray[-50:]
            value = mean(sortarray)
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

def bin2d(nbins,xlist,ylist,zlist,createaxis=True,xaxis=0,yaxis=0):

    if createaxis:
        xmin = min(xlist)
        xmax = max(xlist)
        x_range = xmax-xmin
        xstep = x_range/nbins
        ymin = min(ylist)
        ymax = max(ylist)
        y_range = ymax-ymin
        ystep = y_range/nbins

        xaxis = arange(xmin,xmax+xstep/2,xstep)
        yaxis = arange(ymin,ymax+ystep/2,ystep)

    else:
        xmin = min(xaxis)
        xmax = max(xaxis)
        x_range = xmax-xmin
        xstep = x_range/nbins
        ymin = min(yaxis)
        ymax = max(yaxis)
        y_range = ymax-ymin
        ystep = y_range/nbins


    binarray = []
    for i in range(len(xaxis)):
        i_array = []
        for j in range(len(yaxis)):
            i_array.append([])
        binarray.append(i_array)

    for i in range(len(zlist)):
        xbin = int((xlist[i]-xmin)/xstep)
        ybin = int((ylist[i]-ymin)/ystep)

        if xbin < len(binarray) and ybin < binarray[0] and xbin >= 0 and ybin >= 0:

            binarray[xbin][ybin].append(zlist[i])

    zarray = []
    for i in range(len(binarray)):
        i_array = []
        for j in range(len(binarray[i])):
            if len(binarray[i][j]) > 20:
                sortarray = sort(binarray[i][j])
                sortarray = sortarray[-20:]
                value = mean(sortarray)
                if functions.isnan(value):
                    i_array.append(value)
                else:
                    i_array.append(0)
                
                #i_array.append(max(binarray[i][j]))
                #i_array.append(mean(binarray[i][j]))

            else:
                i_array.append(0)
        zarray.append(i_array)

    zarray = array(zarray)
    zarray = zarray / zarray.max()

    xaxis = xaxis + 0.5*xstep
    yaxis = yaxis + 0.5*ystep
    return zarray,xaxis,yaxis


### Copy the output files
def plot_folder(ax1,ax2,pltcolor):

    folder_params = loadtxt("mcmc_tested_params")
    folder_params = transpose(folder_params)
    folder_prob = loadtxt("mcmc_chisq_log")
    folder_prob = transpose(folder_prob)
    folder_prob = exp((min(folder_prob)-folder_prob)/2)

    #best_params = loadtxt("best_stellar_param_mcmc")

    nlen = min([len(folder_params[0]),len(folder_prob)])

    x = folder_params[ax1][:nlen]
    y = folder_params[ax2][:nlen]
    z = folder_prob[:nlen]

    binarray,xaxis,yaxis = bin2d(30,x,y,z)

    levels = [0.607,0.135]
    #plt.hexbin(y,x,bins=25,cmap="binary",vmax=5)
    plt.scatter(y,x,s=0.1)
    plt.contour(yaxis,xaxis,binarray,levels,colors=pltcolor)

    #plt.scatter(best_params[ax2],best_params[ax1],s=100,color=pltcolor,marker="+")

plt.figure(figsize=(8,8))

param_names = functions.read_table(functions.read_ascii("best_param_mcmc"))[0]

def get_index(pname):
    for i in range(len(param_names)):
        if param_names[i] == pname:
            break
    return i

plt.figure(figsize=(12,8))
### Plot rratio vs rsum
plt.subplot(221)
plot_folder(get_index("rratio"),get_index("rsum"),"r")
plt.xlabel("rsum")
plt.ylabel("rratio")
#plt.show()

### Plot rratio vs i0
plt.subplot(222)
plot_folder(get_index("i_0"),get_index("rratio"),"r")
plt.xlabel("rratio")
plt.ylabel("i0")
#plt.show()

### Plot planet_f planet_alpha
plt.subplot(223)
plot_folder(get_index("planet_f"),get_index("planet_alpha"),"r")
plt.xlabel("planet_alpha")
plt.ylabel("planet_f")
#plt.show()

### Plot planet_f rratio
plt.subplot(224)
plot_folder(get_index("rratio"),get_index("planet_f"),"r")
plt.xlabel("planet_f")
plt.ylabel("rratio")
plt.show()
