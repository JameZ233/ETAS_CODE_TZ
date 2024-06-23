# ETAS_CODE_TZ
 Code for ETAS model

 How the Code works:
 Main File: ETASSimulation.py
 First, setup the parameters. 

 # Set the parameters for the ETAS model
mu = 3.29                             #   Minimum magnitude
K = 1.0                                 #   Aftershock productivity parameter >>>>>>>>>>> Usually have K = 1
alpha = 1.0                             #   Exponent in the productivity relation, multiplies bval
qval = 1.5                              #   Spatial Omori exponent
sigma = 0.5                             #   Could be used in the spatial dependence part (not used in this code)
ntmax =20500                            #   Number of earthquakes requested
mag_large = 7.0                         #   Magnitude of large earthquakes for re-clustering or for stacked aftershock plots
rate_bg = 5.0                           #   Background rate
bval = 0.95                             #   GR b value
pval = 1.2                              #   Omori p-value       >>>>>>>>> Usually pval = 1.2 or so
corr_length = 100                       #   Parameter in the spatial Omori relation (km)
corr_time   = 1                         #   Parameter in the temporal Omori relation (natural time)
dt_ratio_exp = 1.5                      #   Used in the equation that computes time to next event

t = 0.0                                 #   Initial time
kk=0                                    #   Counter for large earthquakes
time_main = 0.                          #   Used in the re-clustering
m  = mu                                 #   Initial magnitude is the minimum magnitude

scale_factor = 0.004                    #   Controls the (re)clustering, this is basically 2 x standard deviation
step_factor_aftrshk = 0.04              #   In degrees, controls the lat-lng steps for the random walk aftershocks

BathNumber = 1.2                        #   Bath's law value.  Note that this can be used to determine the
                                        #       the effective number of aftershocks for a given mainshock magnitude
                                        #       and thus the ratio of aftershocks to mainshocks

scale_mag = 1.                          #   Scaling the aftershocks is from mu+scale_mag up to mag_large
mag_threshold = mu + BathNumber         #   Only events with mags larger than this can have aftershocks
plot_params =  False                     #   Used to display the params in the boxes on the plots
plot_USGS_data  = False                 #   If we want to plot real data, first converts the
                                        #       USGS files to the correct format and file name then runs the functions
                                        
pfit_lo = 80                            #   Low value of the parameters to fit a line to the Omori scaling plot
pfit_hi = 600                           #   High value of the parameters to fit a line to the Omori scaling plot
NSteps  = 36                            #   Used for the Exponential Moving Average in the number-timeseries

#############################################
#############################################

Next, setup the flags for running code and generating output.

generate_and_write_new_catalog              = True   #Creates a new random catalog
write_new_adjusted_file                     = True   #Creates scale invariant aftershocks

plot_magnitude_time                         = False  #Plots the raw catalog
plot_magnitude_time_clustered               = False  #Plots the clustered catalog

plot_magnitude_frequency                    = False  #GR diagram

plot_scale_invariant_omori_plot             = False  #Omori Plot
plot_number_events_vs_time                  = False  #Monthly number of events vs. time for scale invariant catalog

plot_LA_seismicity_map                      = True  #Spatial locations of earthquakes

download_USGS_catalog                       = False  #Input data from USGS

###############################################################################################
###############################################################################################

Last, wait for the plots and simulations with exact locations in USGS_regional.catalog.

Explanation of major functions:
ETASPlotV3.py       #Files that include all the function to make the plots
ETASFileWriter.py   #Files that include the function to write the data saved in text file to csv file
ETASCalcV5.py       #Files that includes functions to do various calculations such as calculating earthquake rate function and cluster function for magnitudes
ETASFileMethods.py   #Files with function that downloads catalog from USGS site

