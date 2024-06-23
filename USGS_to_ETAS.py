import numpy as np
import random
import matplotlib.pyplot as plt
import math
from math import log10
import datetime
import time

import ETASFileWriter

   ###############################################################################################
   
   #    This code writes the USGS catalog into the ETAS catalog format   

   ###############################################################################################
   

def mainshock_cycles(mags, events, label_events, mag_large):

    working_list        =   []     #   Save the first event
    
    #   First count the number of large events in the time series
    
    number_mag_large = 0
    
    for i in range(len(events)):
        if mags[i] >= mag_large:
            number_mag_large += 1
            
    print('number_mag_large: ', number_mag_large)
        
    cycle_list          = [[] for i in range(number_mag_large)]
    label_cycle_list    = [[] for i in range(number_mag_large)]
    
    kk = -1
    for i in range(len(events)):
        if mags[i] >= mag_large:
            kk += 1    
            
#         print('i,kk, number_mag_large, len(events), mags[0]', i, kk, number_mag_large, len(events), mags[0])
        cycle_list[kk].append(events[i])
        label_cycle_list[kk].append(label_events[i])
        
    return cycle_list, label_cycle_list
    
    ###############################################################################################

def omori_aftershocks(cycle_list, events, mags):

    aftershock_list     =   []
    events_adjusted     =   []
    
    #   Set the initial version of aftershock_list = cycle_list
    #       (aftershock_list should be a list of lists where the first 
    #       entry in each list is the time of the mainshock)
    
    aftershock_list = cycle_list
    
    #   Now alter aftershock list so that the intial entry into a list is the time of mainshock,
    #   and the remainder are difference times from the mainshock
    
    for i in range(len(aftershock_list)):
        for j in range(1,len(aftershock_list[i])):      #  First entry is time of mainshock
            aftershock_list[i][j] = aftershock_list[i][j] - aftershock_list[i][0] #   Other entries are difference times
            
    return aftershock_list
    
    ###############################################################################################
    
def write_output_files(date_file, time_file, year_file, x_events, y_events, mags, z_events, aftshk_bg_list,\
        text_file_name, csv_file_name):

    nowcast_output=open(text_file_name, 'w')

    for i in range(len(year_file)):
        print(date_file[i], time_file[i], year_file[i], round(x_events[i],4),round(y_events[i],4),mags[i], \
                round(z_events[i],4), aftshk_bg_list[i], file = nowcast_output)

    nowcast_output.close()
    
    input_file = text_file_name

    output_file = csv_file_name
    
    ETASFileWriter.text_to_csv(input_file, output_file)
    
    return
    
    ###############################################################################################
    
def read_existing_catalog(input_file_name):

    mags                    =   []
    date_events             =   []
    time_events             =   []
    year_events             =   []
    depth_events            =   []
    y_events                =   []
    x_events                =   []

    data_file = open(input_file_name,"r")
    
    for line in data_file:
        items = line.strip().split()
        
        date_events.append(items[0])
        time_events.append(items[1])
        year_events.append(float(items[2]))
        x_events.append(float(items[3]))
        y_events.append(float(items[4]))
        mags.append(float(items[5]))
        depth_events.append(float(items[6]))
        
    data_file.close()  
    
    return date_events, time_events, year_events, x_events, y_events, mags, depth_events
    
    ###############################################################################################
    
def assign_labels_california(counter, nt, m, mu, bval, mag_threshold, BathNumber):
            
    event_occurred = False
        
    if nt == 0:  #   The first event is always larger than the threshold
        counter = int(10** (bval*(m - mu - BathNumber)) * 1.0)
        event_occurred = True
        
    elif m >= mag_threshold and counter <= 0 and nt > 0 and event_occurred == False:    
                                                            #   A large event just occurred that can have aftershocks, we 
                                                            #   so we reset the aftershock counter and pick a 
                                                            #   epicenter at random from the list of possible locations
        counter = int(10** (bval*(m - mu - BathNumber)) * 1.0)
        event_occurred = True
        
    elif m >= mag_threshold and counter > 0 and nt > 0 and event_occurred == False:     
                                                            #   We are in the aftershock train from a large earthquake,
                                                            #   and we just had an aftershocks that can itself have aftershocks
                                                            #   so we we might extend the counter if the aftershock magnitude
                                                            #   was large enough
        new_counter = int(10** (bval*(m - mu - BathNumber)) * 1.0)
        if new_counter > counter:
            counter = new_counter
        event_occurred = True
        
    counter -= 1
    
    return counter
    
    ###############################################################################################
    ###############################################################################################
    


# Set the parameters for the ETAS model
                       #   Background rate - really just a dummy variable
bval = 0.95                             #   GR b value
           #   Used in the re-clustering -- this function not used in this code version
BathNumber = 1.1                        #   From Bath's law

    ###############################################################################################
#       This short code classifies the 
    
input_file_name = './USGS_regional.catalog'
   
date_events, time_events, year_events, x_events, y_events, mags, z_events = \
        read_existing_catalog(input_file_name)
        
counter = 1
mu = min(mags)
mag_threshold = mu + BathNumber         #   Only events with mags larger than this can have aftershocks
mag_large = 7.0

number_background  = 0
number_aftershocks = 0

label_events  =   []

for nt in range(len(mags)):
    m = mags[nt]
    last_counter = counter
    counter = assign_labels_california(counter, nt, m, mu, bval, mag_threshold, BathNumber)
    
    if last_counter < 0:
        label_events.append('b')
        number_background += 1
    else:
        label_events.append('a')
        number_aftershocks += 1
        
print()
fraction_aftershocks = round( 100.0*float(number_aftershocks)/float(number_aftershocks + number_background), 2)
print('Fraction Aftershocks: ', str(fraction_aftershocks) + '%' )

cycle_list, label_cycle_list = mainshock_cycles(mags, year_events,label_events, mag_large)

    
aftershock_list = omori_aftershocks(cycle_list, year_events, mags)

text_file_name = './Output_Files/ETAS_Scale_Invariant_Output.txt'
csv_file_name  = './Output_Files/ETAS_Scale_Invariant_Output.csv'

write_output_files(date_events, time_events, year_events, x_events, y_events, mags, z_events, label_events,\
        text_file_name, csv_file_name)
        
    ###############################################################################################


   
    