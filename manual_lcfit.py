import os
import sys
import functions
import fitting_functions
from numpy import *
import matplotlib.pyplot as plt
from scipy import optimize

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

print "FREE PARAMS"
for i in range(len(free_param_names)):
    print free_param_names[i],free_param_vals[i],free_param_range[i]

print "FIXED PARAMS"
for i in range(len(fixed_param_names)):
    print fixed_param_names[i],fixed_param_vals[i]


fitting_functions.manual_lcfit(free_param_vals,free_param_names,fixed_param_names,fixed_param_vals,lc[0],cadence[0])
