README

Lightcurve fitting script
George Zhou 17/06/2013
george@mso.anu.edu.au

REQUIRES:
python
numpy
scipy 
matplotlib
emcee http://dan.iel.fm/emcee/
gfortran or equivalent

INSTRUCTIONS

1) Edit config_file, input initial parameters. For any parameter, the error defines the space within which MCMC is allowed to explore. Set error to 0 to fix parameter. 
Rsum = (R* + Rp) / a        Rratio = Rp/R*

Under "Fitting options", you can set the number of walkers, burns, mcmc iterations to run. 50, 100, 50000 (respectively) seem to be nice numbers to use, and runs within 15 minutes.

Do not use RUN_DOWNHILL, this seems to be giving the same output as the initial inputs, and is a waste of time - some troubleshooting is required here.

2) run fit algorithm 
python fit_lc.py

This will produce the outputs
mcmc_chisq_log
mcmc_tested_params

You'll need those later for analysis

3) run error analysis
python find_errors.py

This will plot a histogram for each parameter in the MCMC fit. It will print out the 65% confidence region for each parameter.

4) Plot the fit and residuals
python plot_system.py



