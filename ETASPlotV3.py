    #   ---------------------------------------------------------------------------------------
    
    # Copyright 2023 by John B Rundle, University of California, Davis, CA USA
    # 
    # Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
    # documentation files (the     "Software"), to deal in the Software without restriction, including without 
    # limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, 
    # and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
    # 
    # The above copyright notice and this permission notice shall be included in all copies or suSKLantial portions of the Software.
    # 
    # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE 
    # WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
    # COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
    # ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    #
    # Note that some part of these codes were written by the AI program chatGPT, a product of openAI
    #
    #   ---------------------------------------------------------------------------------------

import numpy as np
import random
import math
from math import log10

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import cartopy.io.img_tiles as cimgt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature, LAND, COASTLINE

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.patches as mpatches
from matplotlib import gridspec
from matplotlib.patches import Rectangle

from matplotlib.offsetbox import AnchoredText
from matplotlib.image import imread
import matplotlib.ticker as mticker

import ETASCalcV5

import Compute_EMA_Timeseries

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

   ###############################################################################################
   
def get_params(params):

    mu                          =   params[0]
    K                           =   params[1]
    pval                        =   params[2]
    qval                        =   params[3]
    sigma                       =   params[4]
    ntmax                       =   params[5]
    mag_large                   =   params[6]
    rate_bg                     =   params[7]
    t                           =   params[8]
    kk                          =   params[9]
    time_main                   =   params[10]
    bval                        =   params[11]
    mag_threshold               =   params[12]
    BathNumber                  =   params[13]
    step_factor_aftrshk         =   params[14]
    corr_length                 =   params[15]
    corr_time                   =   params[16]
    alpha                       =   params[17]
    dt_ratio_exp                =   params[18]
    scale_factor                =   params[19]
    plot_params                 =   params[20]
    pfit_lo                     =   params[21]
    pfit_hi                     =   params[22]
    NSteps                      =   params[23]

    return mu, K, pval, qval, sigma, ntmax, mag_large, rate_bg, t, kk, time_main, bval, mag_threshold, \
        BathNumber, step_factor_aftrshk, corr_length, corr_time, alpha, dt_ratio_exp, scale_factor, \
        plot_params, pfit_lo, pfit_hi, NSteps
    
def plot_mag_timeseries_raw(year_events, mags, params):

    mu, K, pval, qval, sigma, ntmax, mag_large, rate_bg, t, kk, time_main, bval, mag_threshold, \
        BathNumber, step_factor_aftrshk, corr_length, corr_time, alpha, dt_ratio_exp, scale_factor, \
        plot_params, pfit_lo, pfit_hi, NSteps\
        = get_params(params)


    mag_list = []
    year_list = []
    
    
    fig, ax = plt.subplots()
    
#     for i in range(int(0.9*len(time_events)), len(time_events)):
    for i in range(int(len(year_events))):
        mag_list.append(mags[i])
        year_list.append(year_events[i])

    ax.plot(year_list, mag_list, marker = 'o', color ='blue', markersize = 1, lw = 0.0)
    
    textstr =   ' Minimum Magnitude: ' + str(mu) +\
                '\n K: ' + str(K) +\
                '\n Alpha: ' + str(K) +\
                '\n Omori p-value: ' + str(pval)+\
                '\n q-value: ' + str(qval)+\
                '\n b-value: ' + str(bval)+\
                '\n Baths Law Number: ' + str(BathNumber)+\
                '\n Aftershock Step: ' + str(step_factor_aftrshk)+ '$^o$'+\
                '\n Correlation Length: ' + str(corr_length)+ ' $Km$'\
                '\n Correlation Time: ' + str(corr_time) + ' Year' +\
                '\n Aftershock Cluster Factor: ' + str(scale_factor) +\
                '\n Rate Ratio Exp:' + str(dt_ratio_exp)

 
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'gray', alpha=0.75)

    # place a text box at lower left
    if plot_params:
        ax.text(0.015, 0.02, textstr, transform=ax.transAxes, fontsize=5,\
            verticalalignment='bottom', horizontalalignment = 'left', bbox=props, linespacing = 1.8)
    
    SupTitle_text = 'Magnitude-Time Diagram for ' + str(len(mags)) +  ' Unclustered ETAS Earthquakes'
    plt.suptitle(SupTitle_text, fontsize=10)
    
    if plot_params:
        Title_text = 'Hyper-Parameter Values Displayed in Legend'
        plt.title(Title_text, fontsize=8)

    plt.xlabel('Time Index')
    plt.ylabel('ETAS Earthquake Magnitude')
    
    FigureName = './Figures/Mag_Time_Raw.png'
          
    plt.savefig(FigureName, dpi=200)
    
    #plt.show()
    
    plt.close()
    
    return
    
    ###############################################################################################
    
def plot_number_timeseries(year_events, mags, params):

    mu, K, pval, qval, sigma, ntmax, mag_large, rate_bg, t, kk, time_main, bval, mag_threshold, \
        BathNumber, step_factor_aftrshk, corr_length, corr_time, alpha, dt_ratio_exp, scale_factor, \
        plot_params, pfit_lo, pfit_hi, NSteps\
        = get_params(params)
   
    month_interval          = 0.07692   #   We use a time interval of "months" 
                                        #       where 1 month = 4 weeks = "lunar month", 52 weeks/year
                                        
    year_index = 0
    max_year = max(year_events)
    min_year = min(year_events)
    
    number_year_bins = int( (max(year_events) - min(year_events))/month_interval)
    
    year_list               =   np.zeros(number_year_bins+1)
    monthly_number_list     =   np.zeros(number_year_bins+1)
    
    for i in range(len(year_list)):
        year_list[i] = min_year + float(i)*month_interval
    
    for i in range(len(year_events)):
        
        last_year_index = year_index
        year_index = int((year_events[i]-min_year)/month_interval)
        if year_index == last_year_index:
            monthly_number_list[year_index]+= 1
            
    log_monthly_number_list = [math.log10(1. + monthly_number_list[i]) for i in range(len(monthly_number_list))]

    fig, ax = plt.subplots()
    
    #   Compute monthly numbers
    
    min_monthly_list = min(log_monthly_number_list)
    
    for i in range(number_year_bins):
        x_eq = [year_list[i],year_list[i]]
        y_eq = [min_monthly_list,log_monthly_number_list[i]]
        ax.plot(x_eq,y_eq, linestyle='-', lw=0.5, color='c', zorder=1, alpha=0.30)
    
    ax.plot(year_list, log_monthly_number_list, marker = 'o', color ='blue', markersize = 1, lw = 0.0, zorder = 5)
    
    #   ---------------------------------------------------
    
    lambda_min_mult = -1.   
    lambda_max_mult = 2.e10 
    
    data_start_year = year_list[0]
#     NSteps = 36
#     NSteps = 12
    year_array = year_list
    
#     timeseries_EMA_reduced, time_list_reduced, log_number_reduced, min_rate, max_rate = \
#             Compute_EMA_Timeseries.build_EMA_timeseries(year_array, NSteps, data_start_year, lambda_min_mult, lambda_max_mult)
#             
    timeseries_EMA_reduced, time_list_reduced, log_number_reduced, min_rate, max_rate = \
            Compute_EMA_Timeseries.build_EMA_timeseries(monthly_number_list, year_list, NSteps, \
            data_start_year, lambda_min_mult, lambda_max_mult)
            
#     print(timeseries_EMA_reduced)
            
    ax.plot(time_list_reduced, log_number_reduced, '--', color ='red', lw = 0.75, zorder = 10, \
            label='EMA, N='+str(NSteps))
            
    title_text = 'EMA Data'
    leg = ax.legend(loc = 'upper right', title=title_text, fontsize=6, fancybox=True, framealpha=0.5)
    leg.set_title(title_text,prop={'size':8})   #   Set the title text font size

    #   ---------------------------------------------------
    
    textstr =   ' Minimum Magnitude: ' + str(mu) +\
                '\n K: ' + str(K) +\
                '\n Alpha: ' + str(K) +\
                '\n Omori p-value: ' + str(pval)+\
                '\n q-value: ' + str(qval)+\
                '\n b-value: ' + str(bval)+\
                '\n Baths Law Number: ' + str(BathNumber)+\
                '\n Aftershock Step: ' + str(step_factor_aftrshk)+ '$^o$'+\
                '\n Correlation Length: ' + str(corr_length)+ ' $Km$'\
                '\n Correlation Time: ' + str(corr_time) + ' Year' +\
                '\n  Cluster Factor: ' + str(scale_factor) +\
                '\n Rate Ratio Exp:' + str(dt_ratio_exp)
                
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'gray', alpha=0.75)

    # place a text box at lower left
    if plot_params:
        ax.text(0.015, 0.02, textstr, transform=ax.transAxes, fontsize=5,\
            verticalalignment='bottom', horizontalalignment = 'left', bbox=props, linespacing = 1.8)
    
    SupTitle_text = 'Monthly Number for ' + str(len(year_events)) +  ' Clustered ETAS Earthquakes'
    plt.suptitle(SupTitle_text, fontsize=10)
    
    if plot_params:
        Title_text = 'Hyper-Parameter Values Displayed in Legend'
        plt.title(Title_text, fontsize=8)

    plt.xlabel('Time (Year)')
    plt.ylabel('$Log_{10}$ (1 + Monthly Number)')
    
    FigureName = './Figures/Number_Time_Clustered.png'
          
    plt.savefig(FigureName, dpi=200)
    
    #plt.show()
    
    plt.close()
    
    return
    
    ###############################################################################################

def plot_mag_timeseries_clustered(year_events, mags, params):


    mu, K, pval, qval, sigma, ntmax, mag_large, rate_bg, t, kk, time_main, bval, mag_threshold, \
        BathNumber, step_factor_aftrshk, corr_length, corr_time, alpha, dt_ratio_exp, scale_factor, \
        plot_params, pfit_lo, pfit_hi, Nsteps\
        = get_params(params)

    mag_list = []
    year_list = []
    
    
    fig, ax = plt.subplots()
    
#     for i in range(int(0.9*len(time_events)), len(time_events)):
    for i in range(int(len(year_events))):
        mag_list.append(mags[i])
        year_list.append(year_events[i])

    ax.plot(year_list, mag_list, marker = 'o', color ='blue', markersize = 1, lw = 0.0)
    
    textstr =   ' Minimum Magnitude: ' + str(mu) +\
                '\n K: ' + str(K) +\
                '\n Alpha: ' + str(K) +\
                '\n Omori p-value: ' + str(pval)+\
                '\n q-value: ' + str(qval)+\
                '\n b-value: ' + str(bval)+\
                '\n Baths Law Number: ' + str(BathNumber)+\
                '\n Aftershock Step: ' + str(step_factor_aftrshk)+ '$^o$'+\
                '\n Correlation Length: ' + str(corr_length)+ ' $Km$'\
                '\n Correlation Time: ' + str(corr_time) + ' Year' +\
                '\n Cluster Factor: ' + str(scale_factor) +\
                '\n Rate Ratio Exp:' + str(dt_ratio_exp)

 
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'gray', alpha=0.75)

    # place a text box at lower left
    if plot_params:
        ax.text(0.015, 0.02, textstr, transform=ax.transAxes, fontsize=5,\
            verticalalignment='bottom', horizontalalignment = 'left', bbox=props, linespacing = 1.8)
    
    SupTitle_text = 'Magnitude-Time Diagram for ' + str(len(mags)) +  ' Clustered ETAS Earthquakes'
    plt.suptitle(SupTitle_text, fontsize=10)
    
    if plot_params:
        Title_text = 'Hyper-Parameter Values Displayed in Legend'
        plt.title(Title_text, fontsize=8)

    plt.xlabel('Time Index')
    plt.ylabel('ETAS Earthquake Magnitude')
    
    FigureName = './Figures/Mag_Time_Clustered.png'
          
    plt.savefig(FigureName, dpi=200)
    
    #plt.show()
    
    plt.close()
    
    return
    
    ###############################################################################################
    
def plot_GR(mags, params):

    mu, K, pval, qval, sigma, ntmax, mag_large, rate_bg, t, kk, time_main, bval, mag_threshold, \
        BathNumber, step_factor_aftrshk, corr_length, corr_time, alpha, dt_ratio_exp, scale_factor, \
        plot_params, pfit_lo, pfit_hi, NSteps\
        = get_params(params)

    mag_bins, freq_mag_bins_pdf, freq_mag_bins_sdf, log_freq_mag_bins = ETASCalcV5.freq_mag(mags)
    freq_mag_bins = [10**(log_freq_mag_bins[i]) for i in range(len(log_freq_mag_bins))]
    
    fig, ax = plt.subplots()

    ax.plot(mag_bins, freq_mag_bins, marker = 'o', color ='blue', markersize = 4, lw = 0)
    
    #   Compute b-value
    mag_lo = mu                 #   log10(10)
    mag_hi = mu+mag_large       #   log10(100)

    mag_line_data, freq_line_data, slope = ETASCalcV5.calc_b_value(log_freq_mag_bins, mag_bins, mag_lo, mag_hi)
    
    bfit = - slope

#   Plot the best fitting line
    ax.plot(mag_line_data, freq_line_data, '--', color = 'red', markersize = 0, lw = 1.5)
    
    ax.set_yscale('log')
    
    textstr1 =   ' Minimum Magnitude: ' + str(mu) +\
                '\n K: ' + str(K) +\
                '\n Alpha: ' + str(K) +\
                '\n Omori p-value: ' + str(pval)+\
                '\n q-value: ' + str(qval)+\
                '\n b-value: ' + str(bval)+\
                '\n Baths Law Number: ' + str(BathNumber)+\
                '\n Aftershock Step: ' + str(step_factor_aftrshk)+ '$^o$'+\
                '\n Correlation Length: ' + str(corr_length)+ ' $Km$'\
                '\n Correlation Time: ' + str(corr_time) + ' Year' +\
                '\n Cluster Factor: ' + str(scale_factor) +\
                '\n Rate Ratio Exp:' + str(dt_ratio_exp)

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'gray', alpha=0.75)

    # place a text box at lower left
    if plot_params:
        ax.text(0.015, 0.02, textstr1, transform=ax.transAxes, fontsize=5,\
            verticalalignment='bottom', horizontalalignment = 'left', bbox=props, linespacing = 1.8)
        
    textstr2 =      ' b-value (from fit): ' + str(round(bfit,2))
                    
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'gray', alpha=0.75)
    ax.text(0.85, 0.975, textstr2, transform=ax.transAxes, fontsize=8,\
        verticalalignment='top', horizontalalignment = 'center', bbox=props, linespacing = 1.8)
 
    SupTitle_text = 'Magnitude-Frequency Diagram for ' + str(len(mags)) +  ' ETAS Earthquakes'
    plt.suptitle(SupTitle_text, fontsize=10)
    
    if plot_params:
        Title_text = 'Hyper-Parameter Values Displayed in Legend'
        plt.title(Title_text, fontsize=8)

    plt.xlabel('Magnitude of ETAS Earthquakes')
    plt.ylabel('$Log_{10}$(Number)')
    
    FigureName = './Figures/GR_Mag_Freq.png'
    plt.savefig(FigureName, dpi=200)

    #plt.show()
    
    plt.close()
    
    return
    

# #     ###############################################################################################
    
def plot_omori_aftershock_decay(aftershock_list, label_cycle_list, \
        number_mainshocks, mag_large_omori, dt_avg, scale_factor, params):
        

    mu, K, pval, qval, sigma, ntmax, mag_large, rate_bg, t, kk, time_main, bval, mag_threshold, \
        BathNumber, step_factor_aftrshk, corr_length, corr_time, alpha, dt_ratio_exp, scale_factor, \
        plot_params, pfit_lo, pfit_hi, NSteps\
        = get_params(params)
        
    bin_length = 5                 #   Days
    bin_length_years = bin_length/365
    
    fig, ax = plt.subplots()
    
    time_bins = np.arange(bin_length, 1000, 1)  #   200 is interestng  -- Is 80* time_bins < min time between M6.5 events??
    
    print('bin length in days: ', bin_length)
    
    aftershock_number_bins = np.zeros(len(time_bins))
    
    for i in range(len(aftershock_list)):
    
        aftershock_sum = ETASCalcV5.sum_aftershock_intrvl(aftershock_list[i])
        
        duration = len(aftershock_list[i])
        
        if aftershock_sum > 0.1:
        
            for j in range(1,len(aftershock_list[i])):
            
                try:
                    bin_index = int( (aftershock_list[i][j])/bin_length_years ) + 1
                except:
                    bin_index = 1
                
                if bin_index < duration:
#                 if bin_index < duration and label_cycle_list[i][j] == 'a':  # If you want to include all events,
#                                                                             #   not just those labeled as aftshocks
                
                    try:
                        aftershock_number_bins[bin_index] += 1
                    except:
                        pass

    print('number_mainshocks', number_mainshocks)    
    
    plot_counter = 0
    for i in range(len(time_bins)):
        j = len(time_bins) - i -1
        if aftershock_number_bins[j] > 0.0:     #   With this catalog, we need to stop at about 80 - maybe the minimum time
                                                #       between mag_large events
            plot_counter += 1
    
#     
    time_bins = [(time_bins[k]*bin_length) for k in range(len(time_bins)) ]
    
#     
    time_bins_reduced               =   []
    aftershock_number_bins_reduced  =   []
    
    time_bins_reduced               = time_bins[2:plot_counter]
    aftershock_number_bins_reduced  = aftershock_number_bins[2:plot_counter]
    
    #   Remove bins that have no events in them
    
    time_bins_to_plot               =   []
    aftershock_number_bins_to_plot  =   []
    
    for i in range(len(time_bins_reduced)):
        if aftershock_number_bins_reduced[i] > 0.001:
            time_bins_to_plot.append(time_bins_reduced[i])
            aftershock_number_bins_to_plot.append(aftershock_number_bins_reduced[i])
            
    ax.plot(time_bins_to_plot, aftershock_number_bins_to_plot, marker = 'o', color = 'blue', markersize = 4, lw = 0)
    
    #   Compute p-value
    time_lo = math.log10(pfit_lo)          
    time_hi = math.log10(pfit_hi)        

    number_data = [math.log10(aftershock_number_bins_to_plot[i]) for i in range(len(time_bins_to_plot))]
    time_data   = [math.log10(time_bins_to_plot[i]) for i in range(len(time_bins_to_plot))]

    time_line_data, number_line_data, slope = ETASCalcV5.calc_p_value(number_data, time_data, time_lo, time_hi)
    
    pfit = - slope

    ax.plot(time_line_data, number_line_data, '--', color = 'red', markersize = 0, lw = 1.5)

#   Plot the best fitting line
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    textstr1 =   ' Minimum Magnitude: ' + str(mu) +\
                '\n K: ' + str(K) +\
                '\n Alpha: ' + str(K) +\
                '\n Omori p-value: ' + str(pval)+\
                '\n q-value: ' + str(qval)+\
                '\n b-value: ' + str(bval)+\
                '\n Baths Law Number: ' + str(BathNumber)+\
                '\n Aftershock Step: ' + str(step_factor_aftrshk)+ '$^o$'+\
                '\n Correlation Length: ' + str(corr_length)+ ' $Km$'\
                '\n Correlation Time: ' + str(corr_time) + ' Year' +\
                '\n Cluster Factor: ' + str(scale_factor) +\
                '\n Rate Ratio Exp:' + str(dt_ratio_exp)

 
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'gray', alpha=0.75)

    # place a text box at lower left
    if plot_params:
        ax.text(0.015, 0.02, textstr1, transform=ax.transAxes, fontsize=5,\
            verticalalignment='bottom', horizontalalignment = 'left', bbox=props, linespacing = 1.8)
            
    textstr2 =      ' Bin Length: ' + str(bin_length) + ' Day(s)'+\
                    '\n p-value (from fit): ' + str(round(pfit,2))
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'gray', alpha=0.75)
    ax.text(0.72, 0.975, textstr2, transform=ax.transAxes, fontsize=8,\
        verticalalignment='top', horizontalalignment = 'left', bbox=props, linespacing = 1.8)
 
    SupTitle_text = 'Aftershock Decay for ' + str(number_mainshocks) +  ' Stacked ETAS Mainshocks with M > ' + str(mag_large_omori)
    plt.suptitle(SupTitle_text, fontsize=10)
    
    if plot_params:
        Title_text = 'Hyper-Parameter Values Displayed in Legend'
        plt.title(Title_text, fontsize=8)

    plt.xlabel('Days After Mainshock')
    plt.ylabel('Number of Aftershocks')
    
    FigureName = './Figures/Omori_Aftershock_Decay.png'

    plt.savefig(FigureName, dpi=200)
    
#     plt.show()
    
    plt.close()
    
    return
    
    ###############################################################################################
    
    
def map_seismicity(NELat_local, NELng_local, SWLat_local, SWLng_local, plot_start_year, \
        Location, catalog,mag_large, min_map_mag, mu, forecast_intervals, \
        date_array, time_array, year_array, lng_array, lat_array, mag_array, depth_array, params):
        

    mu, K, pval, qval, sigma, ntmax, mag_large, rate_bg, t, kk, time_main, bval, mag_threshold, \
        BathNumber, step_factor_aftrshk, corr_length, corr_time, alpha, dt_ratio_exp, scale_factor, \
        plot_params, pfit_lo, pfit_hi, NSteps\
        = get_params(params)

    #   Note:  This uses the new Cartopy interface
    #
    #   -----------------------------------------
    #
    #   Define plot map
    
    dateline_crossing = False
    
    #   Define coordinates
    left_long   = SWLng_local
    right_long  = NELng_local
    top_lat     = SWLat_local
    bottom_lat  = NELat_local
    
    delta_deg_lat = (NELat_local  - SWLat_local) * 0.5
    delta_deg_lng = (NELng_local  - SWLng_local) * 0.5
    
    longitude_labels = [left_long, right_long]
    longitude_labels_dateline = [left_long, 180, right_long, 360]   #   If the map crosses the dateline
    
    central_long_value = 0
    if dateline_crossing:
        central_long_value = 180
        
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=central_long_value))
    ax.set_extent([left_long, right_long, bottom_lat, top_lat])

    # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')

    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='face',
                                        facecolor='coral')
                                        

                                        
    ocean_10m_3000 = cfeature.NaturalEarthFeature('physical', 'bathymetry_H_3000', '10m',
#                                         edgecolor='black',
                                        facecolor='#0000FF',
                                        alpha = 0.3)
                                        

                                        
    lakes_10m = cfeature.NaturalEarthFeature('physical', 'lakes', '10m',
#                                        edgecolor='black',
                                        facecolor='blue',
                                        alpha = 0.75)
                                        
    rivers_and_lakes = cfeature.NaturalEarthFeature('physical', 'rivers_lakes_centerlines', '10m',
#                                        edgecolor='aqua',
                                        facecolor='blue',
                                        alpha = 0.75)

    ax.add_feature(ocean_10m_3000)

    ax.add_feature(cfeature.BORDERS, linestyle='-', alpha=.5, linewidth=0.5)
    ax.add_feature(cfeature.LAKES, alpha=0.95)
    ax.add_feature(cfeature.RIVERS, linewidth= 0.5)
    ax.add_feature(cfeature.STATES, edgecolor='gray',linewidth= 0.5)
#     ax.add_feature(states_provinces, edgecolor='gray')
    ax.coastlines(resolution='10m', color='black', linewidth=0.5)
    
#     stamen_terrain = cimgt.StamenTerrain()
    stamen_terrain = cimgt.Stamen('terrain-background')
    #   Zoom level should not be set to higher than about 6
    ax.add_image(stamen_terrain, 6)

    if dateline_crossing == False:
        gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels=True,
                   linewidth=0.25, color='black', alpha=0.5, linestyle='dotted')
                   
    if dateline_crossing == True:
        gl = ax.gridlines(xlocs=longitude_labels_dateline, draw_labels=True,
                   linewidth=1.0, color='white', alpha=0.5, linestyle='--')

    gl.top_labels = False
    gl.right_labels = False
    gl.xlines = True
    gl.ylines = True

    if catalog == 'LosAngeles':
        gl.xlocator = mticker.FixedLocator([-112,-114,-116, -118, -120, -122])
    
    if catalog == 'Tokyo':
        gl.xlocator = mticker.FixedLocator([132,134,136, 138, 140, 142, 144, 146])

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    
    gl.xlabel_style = {'size': 12, 'color': 'black'}
    gl.ylabel_style = {'size': 12, 'color': 'black'}
    
    #   -----------------------------------------
    #   Put california faults on the map
    
    input_file_name = './California_Faults.txt'
    input_file  =   open(input_file_name, 'r')
    
    for line in input_file:
        items = line.strip().split()
        number_points = int(len(items)/2)
        
        for i in range(number_points-1):
            x = [float(items[2*i]),float(items[2*i+2])]
            y = [float(items[2*i+1]), float(items[2*i+3])]
            ax.plot(x,y,'-', color='darkgreen',lw=0.55, zorder=2)
    
    input_file.close()
    #    
    #   -----------------------------------------
    #
    #   Plot the data
    
    for i in range(len(mag_array)): #   First plot them all as black dots
        if float(mag_array[i]) >= min_map_mag and float(year_array[i]) >= plot_start_year:
            ax.plot(float(lng_array[i]), float(lat_array[i]), '.k', ms=0.5, zorder=1)
        
        if float(mag_array[i]) >= 6.0 and float(mag_array[i]) < 6.89999  and float(year_array[i]) >= plot_start_year:
 #            ax.plot(float(lng_array[i]), float(lat_array[i]), 'g*', ms=11, zorder=2)
            ax.plot(float(lng_array[i]), float(lat_array[i]), 'o', mec='b', mfc='None', mew=1.25, \
                ms=6, zorder=2)
            
        if float(mag_array[i]) >= 6.89999 and float(year_array[i]) >= plot_start_year:
#             ax.plot(float(lng_array[i]), float(lat_array[i]), 'y*', ms=15, zorder=2)
            ax.plot(float(lng_array[i]), float(lat_array[i]), 'o', mec='r', mfc='None', mew=1.25,\
                ms=12, zorder=2)
# 
    lower_mag = min_map_mag

    ax.plot(float(lng_array[i])+ 1000., float(lat_array[i])+1000., '.k', ms=1, \
                    zorder=4, label='$5.9 > M \geq $'+str(lower_mag) )
                    
    ax.plot(float(lng_array[i])+1000., float(lat_array[i])+1000., 'o', mec='b', mfc='None', mew=1.25, \
                    ms=6, zorder=4, label='\n$6.9 > M \geq 6.0$')
                    
    ax.plot(float(lng_array[i])+1000., float(lat_array[i])+1000., 'o', mec='r', mfc='None', mew=1.25,\
                    ms=12, zorder=4, label='\n$M \geq 6.9$')
                    
    title_text = 'Magnitude Key'
    leg = ax.legend(loc = 'lower left', title=title_text, fontsize=6, fancybox=True, framealpha=0.5)
    leg.set_title(title_text,prop={'size':8})   #   Set the title text font size
    
    textstr =   ' Minimum Catalog Magnitude: ' + str(mu) +\
                '\n Minimum Map Magnitude: ' + str(min_map_mag) +\
                '\n K: ' + str(K) +\
                '\n Alpha: ' + str(K) +\
                '\n Omori p-value: ' + str(pval)+\
                '\n q-value: ' + str(qval)+\
                '\n b-value: ' + str(bval)+\
                '\n Baths Law Number: ' + str(BathNumber)+\
                '\n Aftershock Step: ' + str(step_factor_aftrshk)+ '$^o$'+\
                '\n Correlation Length: ' + str(corr_length)+ ' $Km$'\
                '\n Correlation Time: ' + str(corr_time) + ' Year' +\
                '\n Cluster Factor: ' + str(scale_factor) +\
                '\n Rate Ratio Exp:' + str(dt_ratio_exp)

 
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='white', edgecolor = 'gray', alpha=0.75)

    # place a text box at upper right
    if plot_params:
        ax.text(0.985, 0.98, textstr, transform=ax.transAxes, fontsize=5,\
            verticalalignment='top', horizontalalignment = 'right', bbox=props, linespacing = 1.8)

    SupTitle_text = 'Seismicity for $M\geq$' + str(min_map_mag)  + ' after ' + str(plot_start_year)
    plt.suptitle(SupTitle_text, fontsize=14, y=0.98)
#     
    if plot_params:
        Title_text = 'Within ' + str(delta_deg_lat) + '$^o$ Latitude and ' + str(delta_deg_lng) + '$^o$ Longitude of ' + Location +\
            ', Hyper-Parameter Values Displayed in Legend'
        plt.title(Title_text, fontsize=6)
    else:
        Title_text = 'Within ' + str(delta_deg_lat) + '$^o$ Latitude and ' + str(delta_deg_lng) + '$^o$ Longitude of ' + Location 
        plt.title(Title_text, fontsize=8)
    
    
    #   -------------------------------------------------------------

    figure_name = './Figures/Seismicity_Map_' + Location + '_' + str(plot_start_year) + '.png'
    plt.savefig(figure_name,dpi=600)

#     plt.show()
#     
    plt.close()

    #   -------------------------------------------------------------
    
    return None

    ######################################################################