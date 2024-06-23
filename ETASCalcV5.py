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
import matplotlib.pyplot as plt
import math
from math import log10
import datetime
import time

import ETASUtilities
import Compute_EMA_Timeseries


    ###############################################################################################
    
#   Define ETAS functions

 # Define the earthquake productivity function
def A(m, K, bval, alpha):
#     return K*10.0**(0.5*(m-mu))
#     return K*10.0**(0.5*m)
    return K*10.0**(alpha*bval*m)
    
# ...................................................................

# Define the seismicity rate density function
def omega(m):   #   See paper by Abdollah Jalilian  and Jiancang Zhuan,
                #   This is primarily used for the spatial dependence part

    return 1./(sigma*np.sqrt(2*np.pi))*np.exp(-(m-mu)**2/(2*sigma**2))

# ...................................................................


# Define the magnitude-frequency relation
def mag_freq(m, mu,bval):
#     return 10.0**(-1.2*m + 5.5)
    return 10.0**(-bval*(m-mu) + 6.5)

# ...................................................................


# Define the time Omori function
def time_omori(t,catalog, pval, corr_time):

#     corr_time = 10.

#     return ((t-catalog[:,0] + 100.) [catalog[:,0]<t] )**(-pval) 

    return ((t-catalog[:,0] + corr_time) [catalog[:,0]<t] )**(-pval) 

# ...................................................................
    
# Define the space Omori function
def space_omori(t,x,y,catalog, qval, corr_length):

    xpart = ((x - catalog[:,2])**2) [catalog[:,0]<t]
    ypart = ((y - catalog[:,3])**2) [catalog[:,0]<t]
    
#     corr_length = 10.  #   We assume that events are correlated within this distance; this may be the jump distance
    
    try:
        space_omori_term = (math.sqrt(xpart + ypart) + corr_length)**(-2*qval)  #   Constant is to prevent singularity.  
    except:
        space_omori_term = 1.    
        
    return space_omori_term

# ...................................................................

# Define the ETAS earthquake rate function
def lambda_ETAS(m, t, x, y, pval, qval, rate_bg, K, alpha, catalog, catalog_mu, bval, corr_length, corr_time):

    omori_product = time_omori(t,catalog,pval, corr_time)
        
    omori_product_mu = time_omori(t,catalog_mu,pval, corr_time)
    
    space_product = space_omori(t,x,y,catalog, qval, corr_length)
    
    space_product_mu = space_omori(t,x,y,catalog_mu, qval, corr_length)
    
    Term_m  = np.sum( A(catalog[:,1][catalog[:,0]<t], K, bval, alpha) * omori_product * space_product)
    
    Term_mu  = np.sum( A(catalog_mu[:,1][catalog_mu[:,0]<t], K, bval, alpha) * omori_product_mu * space_product_mu)

#     Term_m  = np.sum( A(catalog[:,1][catalog[:,0]<t], K) * omori_product)
#     
#     Term_mu  = np.sum( A(catalog_mu[:,1][catalog_mu[:,0]<t], K) * omori_product_mu)


    rate    =   rate_bg + Term_m
    rate0   =   rate_bg + Term_mu
    
    return rate, rate0
    
# ...................................................................
    
def generate_catalog(params):

#   DefineaArrays and lists
    
    catalog         = np.empty((0, 4))
    catalog_mu      = np.empty((0, 4))

    mags                    = []
    time_events             = []
    rate_vector             = []
    rate_mu_vector          = []
    x_events                = []
    y_events                = []
    z_events                = []

    working_list            = []
    aftershock_list         = []

    Sum_m   =   []
    Sum_mu  =   []
    nt_list =   []
    
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
    
 
    ###############################################################################################
    ###############################################################################################
    
    #   Main program

    #   Initialize variables

    x=0.
    y=0.
    z=10.

    zmin = 0.
    zmax = 20.

    nt = 0
    
    m = mu
    
    step_factor = 1.0
    
    dt_prefactor = 3.5 #   This can be changed to better approximate calendar time, but it 
                        #       is random and unreliable

#   ----------------------------------------------------------------

    #   Read in USGS file for background events
    
    LA_lat = 34.0522
    LA_lng = -188.2437

    USGS_input_file = './USGS_regional.catalog'
    input_file = open(USGS_input_file, 'r')
    
    lat_events  =   []
    lng_events  =   []
    
    lat_lng_all             =   []
    lat_lng_big             =   []
    lat_lng_mid             =   []
    lat_lng_all_Bath        =   []
    
    number_background  = 0
    number_aftershocks = 0
    
    for line in input_file:
        items = line.strip().split()
        mag        = float(items[5])
        if mag >= mu + BathNumber:
            eq_lat = float(items[4])
            eq_lng = float(items[3])
            lat_lng_all_Bath.append((eq_lat,eq_lng))
            
        if mag >= mu:
            eq_lat = float(items[4])
            eq_lng = float(items[3])
            lat_lng_all.append((eq_lat,eq_lng))
            
        if mag >= 5.5 and mag < 6.5:
            eq_lat = float(items[4])
            eq_lng = float(items[3])
            lat_lng_mid.append((eq_lat,eq_lng))
            
        if mag >= 6.5:
            eq_lat = float(items[4])
            eq_lng = float(items[3])
            lat_lng_big.append((eq_lat,eq_lng))
            
    lat_lng_all_Bath = lat_lng_all       
#   ----------------------------------------------------------------

    counter = 1
    
    aftshk_bg_list  =   []
    
    epicenter = random.choice(lat_lng_all)           #   Just pick an epicenter from the existing catalog of California sites
    last_lat = epicenter[0]
    last_lng = epicenter[1]
    
#     data_file = open('data_file_random', 'w')
#     data_file.close()
    


    while nt < ntmax:
    
        nt+= 1
        
#   ----------------------------------------------------------------
    
        length_scale = 1.0
    
#   Choose params for dt, and the spatial random walk

        dt_avg = 1.
    
        if nt >= 500:
            dt_avg = t/float(nt)
            step_factor = dt/dt_avg
        
        rate, rate0 = lambda_ETAS(mu, t, x, y, pval, qval, rate_bg, K, alpha, catalog, catalog_mu, bval, corr_length, corr_time)
    
        ratio_inv = (rate0/rate)**dt_ratio_exp
        
#         dt = - 0.1 * (ratio_inv) *  math.log(1 - ratio_inv * np.random.random()) #  >> try raising ratio_inv**1.5 maybe?
        dt = - dt_prefactor * (ratio_inv) *  math.log(1 - ratio_inv * np.random.random()) #  >> try raising ratio_inv**1.5 maybe?
        
        t += dt
    
        print('nt, t, dt, x, y, m, ratio_inv: ', nt, round(t,4), round(dt,4), round(x,4), round(y,4), round(m,2), round(ratio_inv,3))
        print()
    
        catalog         = np.vstack((catalog, [t, m, x, y]))
        catalog_mu      = np.vstack((catalog_mu, [t, mu, x, y]))  
    
        m = mu + np.random.choice\
                (np.arange(0, 999, 0.01), \
                p=mag_freq(np.arange(0, 999, 0.01), mu, bval)/np.sum(mag_freq(np.arange(0, 990, 0.01), mu, bval)) )
                
        if m > 8.2:
        #   The magnitude must saturate at about m = 8
            m =    8.2*(2.*(  1/(1 + math.exp(-m/2)) - 0.5  ) )
                
                
#         if nt == 500:
#             m = mag_large       #   Always start the file off with a large event
    
        print('event number, t, m, ratio_inv: ', nt, round(t,4),round(m,2), round(ratio_inv,3))
        print('Number of M>=' + str(mag_large) +  ' Events is currently ' + str(kk) )
        print()
        
#   ----------------------------------------------------------------
        
        last_counter = counter
        
        lat_event, lng_event, counter =  \
                lat_lng_all_california(lat_lng_all, lat_lng_all_Bath, lat_lng_mid, lat_lng_big,\
                counter, nt, last_lat, last_lng, m, mu, bval,\
                mag_threshold, step_factor_aftrshk, BathNumber)
                
        last_lat = lat_event
        last_lng = lng_event
        
        #   Convert lng-lat to x,y coordinates relative to LA
        
        x = compute_great_circle_distance(LA_lat, lng_event, LA_lat, LA_lng)
        x = abs(x)
        if lng_event < LA_lng:
            x = -x
        y = compute_great_circle_distance(lat_event, LA_lng, LA_lat, LA_lng)
        y = abs(y)
        if lat_event < LA_lat:
            y = -y

        zrandom = np.random.uniform(-1,1)
        z = z + step_factor*zrandom 

        if z < zmin:
            z = zmin
        elif z > zmax:
            z = zmax
        
#   ----------------------------------------------------------------

        if nt<500:
            t0 = t
    
        if nt >= 500:           #   Neglect any transients
        
            mags.append(m)
            time_events.append(t-t0)
            nt_list.append(nt-500)
    
            rate_vector.append(rate)
            rate_mu_vector.append(rate0)
    
            x_events.append(lng_event)
            y_events.append(lat_event)
            z_events.append(z)
            if last_counter < 0:
                aftshk_bg_list.append('b')
                number_background += 1
            else:
                aftshk_bg_list.append('a')
                number_aftershocks += 1
            
        
#             Sum_m.append(np.sum(A(catalog[:,1][catalog[:,0]<t])* time_omori(t,catalog, pval)))
#     
#             Sum_mu.append(np.sum(A(catalog_mu[:,1][catalog_mu[:,0]<t])* time_omori(t,catalog_mu, pval)))

        print('m, last_lat, last_lng, counter: ', m, last_lat, last_lng, counter)
    
        if m >= mag_large:    #   Write some stuff and accumulate aftershocks for plotting
            
            kk+= 1
            print('******* Event number: ' + str(nt) + ' of ' + str(ntmax) + ' Events at time '+  str(round(t,3)) +\
                ', for event ' + str(kk) + ' that had magnitude M' + str(round(m,3)) + ', which is >= M' + str(mag_large) )
            print()
            time.sleep(1) # Sleep for 2 seconds
            
    print()
    fraction_aftershocks = round( 100.0*float(number_aftershocks)/float(number_aftershocks + number_background), 2)
    print('Fraction Aftershocks: ', str(fraction_aftershocks) + '%' )
    
    return t, dt_avg, mags, time_events, x_events, y_events, z_events, aftshk_bg_list
    
   ###############################################################################################

def mainshock_cycles(mags, events, label_events, mag_large):

    working_list        =   []     #   Save the first event
    
    mags[0] = mag_large #   Just to be sure!
    
    #   First count the number of large events in the time series
    
    number_mag_large = 0
    
    for i in range(len(events)):
        if mags[i] >= mag_large:
            number_mag_large += 1
        
    cycle_list          = [[] for i in range(number_mag_large)]
    label_cycle_list    = [[] for i in range(number_mag_large)]
    
    kk = -1
    for i in range(len(events)):
        if mags[i] >= mag_large:
            kk += 1    
            
        cycle_list[kk].append(events[i])
        label_cycle_list[kk].append(label_events[i])
        
    return cycle_list, label_cycle_list
    
    ###############################################################################################

def omori_recluster(cycle_list, label_cycle_list, events, mags, scale_factor):

    aftershock_list     =   []
    events_adjusted     =   []
    
    #   Set the initial version of aftershock_list = cycle_list
    #       (aftershock_list should be a list of lists where the first 
    #       entry in each list is the time of the mainshock)
    
    aftershock_list = cycle_list
    
    #   Pick a random scale factor centered on the input scale factor for each of the individual cycles
    
    #   Now alter aftershock list so that the intial entry into a list is the time of mainshock,
    #   and the remainder are difference times from the mainshock
    
    for i in range(len(aftershock_list)):
        for j in range(1,len(aftershock_list[i])):      #  First entry is time of mainshock
            aftershock_list[i][j] = aftershock_list[i][j] - aftershock_list[i][0] #   Other entries are difference times

    for i in range(len(aftershock_list)):
    
        scale_factor_random = 1. - scale_factor * np.random.random()
    
        num_exps = 0
        for k in range(1,len(label_cycle_list[i])-1):    # -2 because we don't apply to mainshock time or last time
            if label_cycle_list[i][k] == 'a':
                num_exps += 1

#         deviation = .001*num_exps
#         num_exps = int(num_exps - deviation*np.random.uniform(0,1))
                
        for j in range(1,num_exps):   #   Skip the first entry which is the mainshock time
            
            mult_factor = scale_factor_random**num_exps
                
            delta_events = aftershock_list[i][j]
            delta_events = delta_events *( mult_factor )
            num_exps -= 1
            
            aftershock_list[i][j] = delta_events
            
    for i in range(len(aftershock_list)):
#     
        for j in range(len(aftershock_list[i])):
            if j == 0:
                time_event = aftershock_list[i][j]
            else:
                time_event = aftershock_list[i][j] + aftershock_list[i][0]
            events_adjusted.append(time_event)
            
    return events_adjusted, aftershock_list
    
    ###############################################################################################
    
#   Generate data for frequency magnitude

def freq_mag(mags):

    bin_diff= 0.1
    number_mag_bins = (float(max(mags)) - float(min(mags)) + 1) * 10

    number_mag_bins     =   int(number_mag_bins)
    range_mag_bins      =   int(number_mag_bins)

    freq_mag_bins_pdf       =   np.zeros(number_mag_bins)
    freq_mag_bins_sdf       =   np.zeros(number_mag_bins)
    freq_mag_pdf_working    =   np.zeros(number_mag_bins)
#     mag_array               =   zeros(number_mag_bins)  
    for i in range(len(mags)):
        bin_number      = int(round(((float(mags[i])-float(min(mags)))/float(bin_diff)),1))
        freq_mag_bins_pdf[bin_number]    +=  1

    for i in range(0,range_mag_bins):
        for j in range(i,range_mag_bins):                           # Loop over all bins to compute the GR cumulative SDF
            freq_mag_bins_sdf[i] += freq_mag_bins_pdf[j]

    number_nonzero_bins=0
    for i in range(0,range_mag_bins):                               # Find the number of nonzero bins
        if freq_mag_bins_sdf[i] > 0:
            number_nonzero_bins+=1

    range_bins = int(number_nonzero_bins)                         # Find number of nonzero bins

    log_freq_mag_bins   =   np.zeros(number_nonzero_bins)              # Define the log freq-mag array

    for i in range(0,range_bins):
        if freq_mag_bins_sdf[i] > 0.0:
            log_freq_mag_bins[i] = -100.0                               # To get rid of it.
            log_freq_mag_bins[i] = math.log10(freq_mag_bins_sdf[i])     # Take the log-base-10 of the survivor function

    mag_bins  =   np.zeros(number_nonzero_bins)
    for i in range(0,range_bins):
      mag_bins[i]=min(mags) + float(i)*bin_diff

    #.................................................................
    
    print_data = False

    if print_data:
        print()
        print(mag_bins)
        print()
        print(freq_mag_bins_pdf)
        print()
        print(freq_mag_bins_sdf)
        print()
        print(log_freq_mag_bins)
        print()
        
    return mag_bins, freq_mag_bins_pdf, freq_mag_bins_sdf, log_freq_mag_bins
    
    ###############################################################################################

def decimal_year_to_date(decimal_year):

#     print('decimal_year', decimal_year)

    year = int(decimal_year)
    days_in_year = (datetime.date(year+1, 1, 1) - datetime.date(year, 1, 1)).days
    day_of_year = int((decimal_year - year) * days_in_year)
    date = datetime.datetime.fromordinal(datetime.date(year, 1, 1).toordinal() + day_of_year - 1)
    time_of_day = datetime.timedelta(seconds=int((decimal_year % 1) * 24 * 60 * 60))
    time_of_day = str(time_of_day) +  '.00'
#   return f"{date.strftime('%Y/%m/%d')} {time_of_day}"

    return date.strftime('%Y/%m/%d'), time_of_day
    
    ###############################################################################################
    
def sum_aftershock_intrvl(event_list):

    working_list = []
    
    for i in range(1, len(event_list)):
        working_list.append(event_list[i])
        
#     print('working_list:', working_list)

    aftershock_sum = np.sum(working_list)
    
    return aftershock_sum


    ###############################################################################################
    
def lat_lng_all_california(lat_lng_all, lat_lng_all_Bath, lat_lng_mid, lat_lng_big, counter, nt, last_lat, last_lng, \
            m, mu, bval, mag_threshold, step_factor_aftrshk, BathNumber):
            
#     data_file = open('data_file_random', 'a')

    event_occurred = False
        
    if nt == 0:  #   The first event is always larger than the threshold
        counter = int(10** (bval*(m - mu - BathNumber)) * 1.0)
        epicenter = random.choice(lat_lng_all)           #   Just pick an epicenter from the existing catalog of California sites
        lat_event = epicenter[0]
        lng_event = epicenter[1]
        choice = 0
        event_occurred = True
        
    if m < mag_threshold and counter <= 0 and nt > 0 and event_occurred == False:
        epicenter = random.choice(lat_lng_all)          #   Just pick an epicenter at random from the existing 
                                                        #   catalog of past California sites
        if m >= 5.5 and m < 6.5:
            epicenter = random.choice(lat_lng_mid)
        if m >= 6.5:
            epicenter = random.choice(lat_lng_big)
        lat_event = epicenter[0]
        lng_event = epicenter[1]
        choice = 1
        event_occurred = True
            
    elif m < mag_threshold and counter > 0 and nt > 0 and event_occurred == False:  
                                                        #   Do the random walk, because this small event is still 
                                                        #   an aftershock within the group following the last large event
        #  choose a new site from the random walk
        lat_aftershock, lng_aftershock = random_walk_aftershocks(last_lng, last_lat, step_factor_aftrshk)
        lat_event = lat_aftershock
        lng_event = lng_aftershock
        choice = 2
        event_occurred = True
        
        
    elif m >= mag_threshold and counter <= 0 and nt > 0 and event_occurred == False:    
                                                            #   A large event just occurred that can have aftershocks, we 
                                                            #   so we reset the aftershock counter and pick a 
                                                            #   epicenter at random from the list of possible locations
        counter = int(10** (bval*(m - mu - BathNumber)) * 1.0)
        
        if m >= mu+BathNumber:
            epicenter = random.choice(lat_lng_all_Bath)     
        if m >= 5.5 and m < 6.5:
            epicenter = random.choice(lat_lng_mid)
        if m >= 6.5:
            epicenter = random.choice(lat_lng_big)
            
        lat_event = epicenter[0]
        lng_event = epicenter[1]
        choice = 3
        event_occurred = True
        
    elif m >= mag_threshold and counter > 0 and nt > 0 and event_occurred == False:     
                                                            #   We are in the aftershock train from a large earthquake,
                                                            #   and we just had an aftershocks that can itself have aftershocks
                                                            #   so we we might extend the counter if the aftershock magnitude
                                                            #   was large enough
        new_counter = int(10** (bval*(m - mu - BathNumber)) * 1.0)
        if new_counter > counter:
            counter = new_counter
        #  choose a new site from the random walk
        lat_aftershock, lng_aftershock = random_walk_aftershocks(last_lng, last_lat, step_factor_aftrshk)
        lat_event = lat_aftershock
        lng_event = lng_aftershock
        choice = 4
        event_occurred = True
        
    last_lng = lng_event
    last_lat = lat_event
    counter -= 1
    


    
#     print(m, lat_event, lng_event, last_lat, last_lng, counter, choice, file = data_file)
#     data_file.close()

    return lat_event, lng_event, counter
    
    ###############################################################################################
    
def random_walk_aftershocks(last_lng, last_lat, step_factor_aftrshk):

    lng_random = np.random.uniform(-1,1)
    lat_random = np.random.uniform(-1,1)

    lng_aftershock = last_lng + step_factor_aftrshk * lng_random
    lat_aftershock = last_lat + step_factor_aftrshk * lat_random      

    return  lat_aftershock, lng_aftershock
    
    ###############################################################################################
    
def compute_great_circle_distance(lat_1, lng_1, lat_2, lng_2):

    # Build an array of x-y values in degrees defining a circle of the required radius

    pic = 3.1415926535/180.0
    Radius = 6371.0

    lat_1 = float(lat_1) * pic
    lng_1 = float(lng_1) * pic
    lat_2 = float(lat_2) * pic
    lng_2 = float(lng_2) * pic
    
    delta_lng = lng_1 - lng_2
    
    delta_radians = math.sin(lat_1)*math.sin(lat_2) + math.cos(lat_1)*math.cos(lat_2)*math.cos(delta_lng)
    if delta_radians > 1.0:
        delta_radians = 1.0
    delta_radians = math.acos(delta_radians)
    
    great_circle_distance = delta_radians * Radius

    return great_circle_distance

    ##################################################################
    
def calc_p_value(number_data, time_data, time_lo, time_hi):
    
    #   Here us where we count the number of data to fit, then define xr[:] and yr[:]
    #   All done in linear corrdinates
    
    print_data = True
    
    xr    =   []
    yr    =   []

    k = -1
    
    for i in range(len(time_data)):
        if (time_data[i] >= time_lo) and (time_data[i] <= time_hi):
            k = k + 1
            xr.append(time_data[i])
            yr.append(number_data[i])

    kmax = k+1
  
    if kmax > 1:
        slope, cept, errs, errc, s2 = ETASUtilities.linfit(xr,yr, kmax)
        if print_data:
            print(' ')
            print('  p-value for the region is:')
            print(-1*round(slope,3),' +/- ',round(errs,3))
            print(' ')
            print(' Intercept for the region is: ', round(cept,3))
            print(' ')
            
    xline = []
    yline = []
    
    for i in range(len(xr)):
        xterm = slope * ( xr[i] - xr[0] ) 
        xline.append(xr[i])  #   slope is negative
        yline.append(yr[0] +  xterm)
        
    #   Now minimize the bias of the line relative to the actual data, i.e., compute a constant
    #       correction to the y-values
    
    error   =   []
    for i in range(len(xr)):
        err_val = (yline[i] - yr[i])
        error.append(err_val)
    mean_error = np.mean(error)
    
#     
#     
    yline   =   [10**(yline[i]- mean_error) for i in range(len(yline))]     
    xline   = [10**(xline[i]) for i in range(len(xline))]

    return xline, yline, slope
            
    ##################################################################
    
def calc_b_value(log_freq_data, mag_data, mag_lo, mag_hi):
    
    #   Here us where we count the number of data to fit, then define xr[:] and yr[:]
    
    print_data = True
    
    xr    =   []
    yr    =   []

    k = -1
    
    for i in range(len(mag_data)):
        if (mag_data[i] >= mag_lo) and (mag_data[i] <= mag_hi):
            k = k + 1
            xr.append(mag_data[i])
            yr.append(log_freq_data[i])

    kmax = k+1
    
  
    if kmax > 1:
        slope, cept, errs, errc, s2 = ETASUtilities.linfit(xr,yr, kmax)
        if print_data:
            print(' ')
            print('  b-value for the region is:')
            print(-1.*round(slope,3),' +/- ',round(errs,3))
            print(' ')
            print(' Intercept for the region is: ', round(cept,3))
            print(' ')
            
    
    xline = []
    yline = []
    
    for i in range(len(xr)):
        xterm = slope * ( xr[i] - xr[0] ) 
        xline.append(xr[i])  #   slope is negative
        yline.append(yr[0] + xterm)
        
    #   Now minimize the bias of the line relative to the actual data, i.e., compute a constant
    #       correction to the y-values, and return the yline value to its expontiated form
    
    error   =   []
    for i in range(len(xr)):
        error.append(yline[i]-yr[i])
    mean_error = np.mean(error)
    
    yline   =   [10**(yline[i]-mean_error) for i in range(len(yline))]        
        
    return xline, yline, slope
            
    ##################################################################



