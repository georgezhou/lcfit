from numpy import *
import matplotlib.pyplot as plt
import functions,sys

def cutout_regions(data,input_regions):
    newdata = []
    for region in input_regions:
        for i in data:
            if i[0] > region[0] and i[0] < region[1]:
                newdata.append(i)
    newdata = array(newdata)
    return newdata
    

data = loadtxt("kplr006603043-2011024051157_slc.fktrap")
param = transpose(loadtxt("param.txt"))

data_midtime = median(data[:,0])
epoch = round((data_midtime-param[0])/ param[1])
midtime = param[0]+epoch*param[1]
tdur = param[2]/24.
oot = [[midtime-tdur*4,midtime-tdur*0.6],[midtime+tdur*0.6,midtime+tdur*4]]
transit = [[midtime-tdur*1.,midtime+tdur*1.]]


# oot = cutout_regions(data,oot)
# plt.scatter(oot[:,0],oot[:,2],s=1)
# oot_fit = polyfit(oot[:,0],oot[:,2],6)
# x = arange(min(oot[:,0]),max(oot[:,0]),0.05)
# oot_plot = polyval(oot_fit,x)
# plt.plot(x,oot_plot)
# plt.show()

data = cutout_regions(data,transit)
# data[:,2] = data[:,2] - polyval(oot_fit,data[:,0])

mag = data[:,2]
flux = 10**((mag-median(mag))/-2.5)

o = open("1.dat","w")
output = [data[:,0],flux,ones(len(data))]
output = transpose(output)
functions.write_table(output,o)
o.close()

plt.scatter(data[:,0],flux)
plt.show()
