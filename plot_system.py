import os
import sys
import functions
import fitting_functions
from numpy import *
import matplotlib.pyplot as plt
from scipy import optimize

def get_maxmin(phase):
    start = []
    end = []
    for i in phase:
        if i > 0.5:
            end.append(i)
        else:
            start.append(i)
    return min(end),1+max(start)


def smooth(x,window_len=11,window='flat'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=ones(window_len,'d')
    else:
        w=eval(window+'(window_len)')

    y=convolve(w/w.sum(),s,mode='valid')

    y = y[window_len/2:-1*window_len/2+1]

    return y

def movingaverage(interval, window_size):
    window= ones(int(window_size))/float(window_size)
    return convolve(interval, window, 'same')

### Global constants
au = 1.496*10**11
msun = 1.988435*10**30
rsun = 6.955*10**8
mjup = 1.8988*10**27 
rjup = 6.9173*10**7
day = 60.*60.*24.
gconst = 6.67*10**(-11)


### Read from config file
### Load initial parameters and lightcurve

temp_param_names = []
temp_param_vals = []
temp_param_range = []

lclist = functions.read_ascii(functions.read_config_file("INPUT_LC_LIST"))
lc = []
for lc_n in lclist:
    lc.append(loadtxt(lc_n))

cadence = []
for i in lc:
    cadence.append(functions.find_cadence(i))

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
temp_param_range.append(eval(functions.read_config_file("PERIOD_ERR")))

temp_param_names.append("t0")
temp_param_vals.append(float(functions.read_config_file("T0"))-floor(float(functions.read_config_file("T0"))))
temp_param_range.append(eval(functions.read_config_file("T0_ERR")))

temp_param_names.append("rsum")
temp_param_vals.append(float(functions.read_config_file("RSUM")))
temp_param_range.append(eval(functions.read_config_file("RSUM_ERR")))

temp_param_names.append("rratio")
temp_param_vals.append(float(functions.read_config_file("RRATIO")))
temp_param_range.append(eval(functions.read_config_file("RRATIO_ERR")))

temp_param_names.append("i_0")
temp_param_vals.append(float(functions.read_config_file("I0")))
temp_param_range.append(eval(functions.read_config_file("I0_ERR")))

temp_param_names.append("ecosw")
temp_param_vals.append(float(functions.read_config_file("ECOSW")))
temp_param_range.append(eval(functions.read_config_file("ECOSW_ERR")))

temp_param_names.append("esinw")
temp_param_vals.append(float(functions.read_config_file("ESINW")))
temp_param_range.append(eval(functions.read_config_file("ESINW_ERR")))

temp_param_names.append("edepth")
temp_param_vals.append(float(functions.read_config_file("EDEPTH")))
temp_param_range.append(eval(functions.read_config_file("EDEPTH_ERR")))

temp_param_names.append("beta")
temp_param_vals.append(float(functions.read_config_file("BETA")))
temp_param_range.append(eval(functions.read_config_file("BETA_ERR")))

temp_param_names.append("fratio")
temp_param_vals.append(float(functions.read_config_file("FRATIO")))
temp_param_range.append(eval(functions.read_config_file("FRATIO_ERR")))

temp_param_names.append("theta")
temp_param_vals.append(float(functions.read_config_file("THETA")))
temp_param_range.append(eval(functions.read_config_file("THETA_ERR")))

temp_param_names.append("phi")
temp_param_vals.append(float(functions.read_config_file("PHI")))
temp_param_range.append(eval(functions.read_config_file("PHI_ERR")))

temp_param_names.append("Protot")
temp_param_vals.append(float(functions.read_config_file("ROTP")))
temp_param_range.append(eval(functions.read_config_file("ROTP_ERR")))

temp_param_names.append("planet_f")
temp_param_vals.append(float(functions.read_config_file("PLANET_F")))
temp_param_range.append(eval(functions.read_config_file("PLANET_F_ERR")))

temp_param_names.append("planet_alpha")
temp_param_vals.append(float(functions.read_config_file("PLANET_ALPHA")))
temp_param_range.append(eval(functions.read_config_file("PLANET_ALPHA_ERR")))

temp_param_names.append("tdiff")
temp_param_vals.append(float(functions.read_config_file("TDIFF")))
temp_param_range.append(eval(functions.read_config_file("TDIFF_ERR")))


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

print "FREE PARAMS"
for i in range(len(free_param_names)):
    print free_param_names[i],free_param_vals[i],free_param_range[i]


print "FIXED PARAMS"
for i in range(len(fixed_param_names)):
    print fixed_param_names[i],fixed_param_vals[i]

x0 = zeros(len(free_param_names))

free_param_vals = [functions.read_ascii("best_param_mcmc")[1]]
free_param_vals = array(functions.read_table(free_param_vals))[0]

print free_param_vals

phase,flux,err,model = [],[],[],[]


for n in range(len(lclist)):


    lc = lclist[n]
    lc = loadtxt(lc)
    phase_n,flux_n,err_n,model_n = fitting_functions.lc_chisq(free_param_vals,free_param_names,fixed_param_names,fixed_param_vals,lc,False,True,cadence[n])
    phase += list(phase_n)
    flux += list(flux_n)
    err += list(err_n)
    model += list(model_n)

####################
### Plot transit ###
####################
phase,flux,err,model = array(phase),array(flux),array(err),array(model)

# ### Test create data
# noise = []
# import random
# for i in phase:
#     noise.append(random.gauss(1,0.00005))
# noise = array(noise)
# #noise = ones(len(phase))

# o = open("temp.txt","w")
# functions.write_table(transpose([phase*110.321612962+1030.3645,model*noise,ones(len(phase))]),o)
# o.close()

# d = loadtxt("temp.txt")
# plt.scatter(phase,d[:,1])
# plt.show()

# sys.exit()

### #### ###

x1,x2 = get_maxmin(phase)

plt.figure(figsize=(8,10))


### Sort in phase
data = transpose([phase,flux,err,model])
data = functions.sort_array(data,0)
phase,flux,err,model = transpose(data)

### Plot data
plt.clf()
plt.subplot(211)
#plt.title("f=0.05, alpha=45")

plt.scatter(phase,flux,s=1,color="k")
plt.scatter(phase+1,flux,s=1,color="k")

plt.scatter(phase,model,s=2,color="r")
plt.scatter(phase+1,model,s=2,color="r")

plt.xlim(x1,x2)
#plt.show()

### Plot residual
#plt.clf()

plt.subplot(212)
residual = flux-model


residual_mvavg = smooth(residual,window_len=15*(len(lclist)),window='flat')
#residual_mvavg = smooth(residual,window_len=60,window='flat')

plt.scatter(phase,residual,s=1,color="k")
plt.scatter(phase+1,residual,s=1,color="k")

#plt.scatter(phase,residual_mvavg,s=10,color="r")
#plt.scatter(phase+1,residual_mvavg,s=10,color="r")
plt.plot(phase,residual_mvavg,"b-")
plt.plot(phase+1,residual_mvavg,"b-")


plt.axhline(y=0,color="r")

plt.xlim(x1,x2)


plt.show()




# ####################
# ### Plot eclipse ###
# ####################

# phase,flux,err,model = array(phase),array(flux),array(err),array(model)


# ### Sort in phase
# data = transpose([phase,flux,err,model])
# data = functions.sort_array(data,0)
# phase,flux,err,model = transpose(data)

# ### Plot data
# plt.clf()
# plt.scatter(phase,flux,s=1,color="k")
# plt.scatter(phase+1,flux,s=1,color="k")

# plt.scatter(phase,model,s=2,color="r")
# plt.scatter(phase+1,model,s=2,color="r")

# plt.xlim(0.55,0.64)
# plt.ylim(0.9995,1.0005)
# plt.show()

# ### Plot residual
# plt.clf()

# residual = flux-model
# #residual_mvavg = movingaverage(residual,60)
# #trend = smooth(residual,window_len=60,window='hamming')
# #residual = residual - trend

# residual_mvavg = smooth(residual,window_len=12,window='flat')
# #residual_mvavg = smooth(residual,window_len=60,window='flat')



# plt.scatter(phase,residual,s=1,color="k")
# plt.scatter(phase+1,residual,s=1,color="k")

# #plt.scatter(phase,residual_mvavg,s=10,color="r")
# #plt.scatter(phase+1,residual_mvavg,s=10,color="r")
# plt.plot(phase,residual_mvavg,"b-")
# plt.plot(phase+1,residual_mvavg,"b-")


# plt.axhline(y=0,color="r")

# plt.xlim(0.55,0.64)
# plt.show()
