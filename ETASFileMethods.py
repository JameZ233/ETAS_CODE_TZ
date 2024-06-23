#!/opt/local/bin python
#import sys
#sys.path.reverse()

    #   Earthquake Methods library of methods and functions
    #   
    #   This code base collects the methods and functions used to make
    #   plots and maps of earthquake data and activity
    #
    ######################################################################

import sys
import matplotlib
import matplotlib.mlab as mlab
from matplotlib import cbook
#from matplotlib.pyplot import figure, show
import numpy as np
from numpy import *
# from mpl_toolkits.basemap import Basemap
from array import array
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.patches as patches
from textwrap import dedent

import scipy.stats
from scipy.stats import poisson

import datetime
import dateutil.parser
from datetime import datetime

import urllib.request, urllib.parse, urllib.error
import urllib.request, urllib.error, urllib.parse
import os

import math
from math import exp

from matplotlib import cm

import http.client
from urllib.error import HTTPError

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

from tabulate import tabulate

import time

    ###########################################################################
    ###########################################################################

    # Get the catalog parameters
    
def get_current_date_time():

    now = datetime.now()

    print("Current year: %d" % now.year)
    print("Current month: %d" % now.month)
    print("Current day: %d" % now.day)
    print("Current hour: %d" % now.hour)
    print("Current minute: %d" % now.minute)
    print("Current second: %d" % now.second)
    
    slash = '/'
    colon = ':'
    
    year  = str(now.year)
    month = str(now.month)
    day   = str(now.day)
    
    if now.month < 10:
        month = '0'+ str(now.month)
        
    if now.day < 10:
        day = '0'+ str(now.day)
    
    current_date = year + slash + month + slash + day
    current_time = str(now.hour) + colon + str(now.minute) + colon + str(now.second)
    
    return (current_date, current_time, year, month, day)
    
def download_base_catalog(NELat, NELng, SWLat, SWLng,Magnitude,begin_date,end_date,output_file_name):

    #   For instructions, refer to:  https://earthquake.usgs.gov/fdsnws/event/1/

    data = {
    
    # Adjust these  -   CA-NV
        "minmagnitude": Magnitude,
        "minlatitude": SWLat,
        "maxlatitude": NELat,
        "minlongitude": SWLng,
        "maxlongitude": NELng,
        "mindepth": -10.0,         # Leaving these depth params in leads to missing some earthquake events
        "maxdepth": 1000.0,


    # Date parameters
#        "starttime": "2070/01/01",
        "starttime": begin_date,
        "endtime": end_date,


    # Leave these
        "eventtype": "earthquake",
        "format": "csv",
        "orderby": "time-asc",
    }
    
    block_size = 10000
    event_offset = 1
#     event_count = [0 for i in range(1)]
    
   
    #   First do a count of events
    
    url = "https://earthquake.usgs.gov/fdsnws/event/1/count?"
    params = urllib.parse.urlencode(data)
    query_string = url + str(params)
    
    error_code = 0
    
    try:
        response_count = urllib.request.urlopen(query_string)
        event_count = response_count.readlines()
#         print('event_count: ', event_count)
        event_count = event_count[0].decode('UTF-8')
#         print('event_count: ', event_count)
        number_events = int(event_count)
        
    except:
        error_code = 1
        number_events = 0
        print('')
        print('1 Download barfed, Error Code: ', error_code)
        pass

    print('')
    print('Number of Events: ', number_events)
    n_blocks = int(event_count)/block_size        #   Number of complete blocks
    n_blocks = int(n_blocks)
    print('Number of complete blocks of size ' + str(block_size) + ' =', n_blocks)
    print('')
    
    for i_block in range(0,n_blocks):
    
    #   -------------------------------------------------------------
    
        event_offset = i_block * block_size + 1
    
        data.update({'offset':event_offset})
        data.update({'limit':block_size})
    
        url = "https://earthquake.usgs.gov/fdsnws/event/1/query?"
        params = urllib.parse.urlencode(data)
        query_string = url + str(params)
    
        error_code = 0
        
        try:
            response = urllib.request.urlopen(query_string)
            catalog = response.readlines()
        except:
            error_code = 1
            print('')
            print('2 Download barfed, Error Code: ', error_code)
            pass
        
        print('Block Number: ', i_block)
        
    #   -------------------------------------------------------------
    
        write_to_file(output_file_name,catalog)
                
    residual_events = number_events
    
    if number_events > block_size:
        residual_events = number_events%(n_blocks*block_size)
        
    if residual_events > 0:
        if n_blocks == 0:
            event_offset = 1
        if n_blocks > 0:
            event_offset = n_blocks * block_size + 1
        data.update({'offset':event_offset})
        data.update({'limit':block_size})
    
        url = "https://earthquake.usgs.gov/fdsnws/event/1/query?"
        params = urllib.parse.urlencode(data)
        query_string = url + str(params)
    
        error_code = 0
        try:
            response = urllib.request.urlopen(query_string)
            catalog = response.readlines()
        except:
            error_code = 1
            print('')
            print('3 Download barfed, Error Code: ', error_code)
            pass
        
        if error_code == 0:
            write_to_file(output_file_name,catalog)
    
    return None
    
    #############################################################
    
def get_base_catalog(NELat, NELng, SWLat, SWLng, completeness_mag, start_date):

    time_seconds = 0.0
    restart = 'NO'          #   We want to output a new file
    output_file_name = "USGS_Base.catalog"
    
    if restart == 'NO':     #   Erase the prior contents of the file
        output_file = open(output_file_name, "w")   #   This statement dumps any previous file
        output_file.close()
        
    current_date, current_time, current_year, current_month, current_day = get_current_date_time()
    print('current_date, current_time: ', current_date, current_time)
    
    start_year = int(start_date.split('/')[0])
    end_year = int(current_year)    
    number_years = end_year - start_year + 1
#
#   -------------------------------------------------------
#     
#   Write the WorldWide catalog by appending data year by year

    print('')
    print('------------------------------------------')
    print('')
    print('Downloading USGS Master Catalog for M>' + str(completeness_mag))
    print('')
    print('------------------------------------------')
    print('')
    
    for i in range (0,number_years):
        begin_year = str(int(i) + start_year)
        begin_date = str(begin_year)+ '/01/01'
        end_date   = str(int(begin_year)+1) + '/01/01'   
        last_date  = str(begin_year) + '/12/31'
        
        print('')
        print('')
        print('Begin Date to End Date: ', begin_date, ' to ', last_date) 

        download_base_catalog(NELat, NELng, SWLat, SWLng, completeness_mag, begin_date,end_date,output_file_name)
        
        time.sleep(time_seconds)
        
    return
    
    #############################################################

def write_to_file(output_file_name,catalog):

    #   This function writes the earthquake data to file
    
#   output_file_name = "USGS_WorldWide.catalog"
    output_file = open(output_file_name, "a")
    
    #   Write output record file

    # Format the data

    date_string = ' '
    time_string = ' '
    ts = 0.
    
    i=-1
    for line in catalog:
        i += 1
        if i > 0:
            line_decode = line.decode('UTF-8')  #   Convert from Byte literal to String
            items = line_decode.strip().split(',')

            date = dateutil.parser.parse(items[0].split('T')[0])     # gives year-month-day (space) hours:minutes:seconds

            ts = date.year + float(date.strftime("%j"))/366

            lat                 = items[1]
            lon                 = items[2]
            dep                 = items[3]
            mag                 = items[4]
            dsb                 = items[0].split('T')[0].split('-')
            date_string         = dsb[0]+'/'+ dsb[1] + '/' + dsb[2]
            time_string         = items[0].split('T')[1]
            time_string         = time_string[:-1]
            
    #   These next checks are in case depths are missing, or depths & mags are listed as 'Unk'

            if dep == 'Unk':
                dep = '0.0'

            if mag == 'Unk':
                mag = '0.0'

            if mag == 'SKL':
                mag = items[4]
                dep = '0.0'

            if mag == 'Mb':
                mag = items[4]
                dep = '0.0'

            if mag == 'Mw':
                mag = items[4]
                dep = '0.0'

            if mag == 'Mc':
                mag = items[4]
                dep = '0.0'

            if mag == 'Md':
                mag = items[4]
                dep = '0.0'

            if mag == 'Ms':
                mag = items[4]
                dep = '0.0'
                
            output_file.write("%s\t%s\t%f\t%s\t%s\t%s\t%s\n" % (date_string, time_string, ts, lon, lat, mag, dep))
                
    output_file.close()
                
    return None

    #################################################################
    
def get_regional_catalog(NELat_local, NELng_local, SWLat_local, SWLng_local, minimum_mag, max_depth,\
        region_catalog_date_start, region_catalog_date_end):
        
#   This code will use shapely.py to filter the base catalog
#       into a rectangular region


    data_string = 'Building catalog for local region...'
    
    print('')
    print('------------------------------------------')
    print('')
    print(data_string)
    print('')
    print('------------------------------------------')
    print('')

    number_polygon_vertices = 4

    #   Construct the string of polygon vertices.  Note that the order is lat, long pairs
    
    vertex_lat = []
    vertex_lng = []
    
    #   Order of vertices of large rectangular region:  NW, NE, SE, SW
    
    vertex_lat.append(NELat_local)
    vertex_lat.append(NELat_local)
    vertex_lat.append(SWLat_local)
    vertex_lat.append(SWLat_local)
    
    vertex_lng.append(SWLng_local)
    vertex_lng.append(NELng_local)
    vertex_lng.append(NELng_local)
    vertex_lng.append(SWLng_local)
    
    point_list = []

    for i in range(number_polygon_vertices):
        point_list.append((float(vertex_lat[i]),float(vertex_lng[i])))
    
    polygon = Polygon(point_list)
    
#
#   -------------------------------------------------------
#    
    
    input_file_name     = "USGS_Base.catalog"
    input_file          =  open(input_file_name, "r")

    output_file_name = "USGS_regional.catalog"
    output_file = open(output_file_name, "w")
    
    for line in input_file:
        items = line.strip().split()
        dep    = items[6]
        mag    = items[5]
        eq_lat = items[4]
        eq_lng = items[3]
        
        point = Point((float(eq_lat),float(eq_lng)))
        
        if (float(dep) <= float(max_depth) and float(mag) >= float(minimum_mag) and polygon.contains(point) == True):
            print(items[0],items[1],items[2],items[3],items[4],items[5],items[6], file=output_file)
        
    output_file.close()
    input_file.close()

    return
    
    #################################################################
        
def read_regional_catalog(min_mag):

    mag_array   =   []
    date_array  =   []
    time_array  =   []
    year_array  =   []
    depth_array =   []
    lat_array   =   []
    lng_array   =   []

    data_file = open("USGS_regional.catalog","r")
    
    for line in data_file:
        items = line.strip().split()
        
        try:
        
            lat                 = items[4]
            lon                 = items[3]
            dep                 = items[6]
            mag                 = items[5]
            date_string         = items[0]
            time_string         = items[1]
            ts                  = items[2]

            if float(mag) >= float(min_mag):
                mag_array.append(mag)           #   List of magnitudes
                date_array.append(date_string)
                time_array.append(time_string)
                year_array.append(ts)
                depth_array.append(dep)
                lat_array.append(lat)
                lng_array.append(lon)
            
        except:
            pass
    
    data_file.close()  

    return mag_array, date_array, time_array, year_array, depth_array, lat_array, lng_array

    ######################################################################
    
   
def read_grid_file():

    #   Short code to read eigenvector file.  Returns eigenvectors in descending order
    #       of eigenvalues

    count = 0
    input_file = open('gridboxes.txt','r')
    for line in input_file:
        count += 1
    input_file.close()
    
    lat         =   [] 
    lng         =   []
    lat_index   =   []
    lng_index   =   []
    
    input_file = open('gridboxes.txt','r')
    for line in input_file:
        items = line.strip().split()
        lat.append(float(items[0]))
        lng.append(float(items[1]))
        lat_index.append(int(items[2]))
        lng_index.append(int(items[3]))
        
    input_file.close()            

    return lat, lng, lat_index, lng_index
    
    ######################################################################
    
def get_timeseries_data(min_mag):

    mag_array, date_array, time_array, year_array, depth_array, lat_array, lng_array = \
            read_regional_catalog(min_mag)

    #   First read the timeseries.txt file.  First determine number of time series
    
    input_file = open('timeseries.txt','r')
    number_ts=0
    for line in input_file:
        number_ts += 1
    input_file.close()
    
    time_bins   = []
    timeseries  = [[] for i in range(number_ts-1)]
    
    input_file = open('timeseries.txt','r')
    i=0
    for line in input_file:
        items = line.strip().split()

        if i == 0:
            for j in range(len(items)):
                time_bins.append(float(items[j]))
        else:
            for j in range(len(items)):
                timeseries[i-1].append(float(items[j]))
                
        i += 1
                
    input_file.close()        

    return time_bins, timeseries
    
    ######################################################################
    
def get_timeseries_EMA_data(min_mag):

    mag_array, date_array, time_array, year_array, depth_array, lat_array, lng_array = \
            read_regional_catalog(min_mag)

    #   First read the timeseries.txt file.  First determine number of time series
    
    input_file = open('timeseries_EMA.txt','r')
    number_ts=0
    for line in input_file:
        number_ts += 1
    input_file.close()
    
    time_bins   = []
    timeseries_EMA = [[] for i in range(number_ts-1)]
    
    input_file = open('timeseries_EMA.txt','r')
    i=0
    for line in input_file:
        items = line.strip().split()

        if i == 0:
            for j in range(len(items)):
                time_bins.append(float(items[j]))
        else:
            for j in range(len(items)):
                timeseries_EMA[i-1].append(float(items[j]))
                
        i += 1
                
    input_file.close()        

    return time_bins, timeseries_EMA
    
    ######################################################################
    
def get_timeseries_EMA_data_LS(LS, min_mag):

    mag_array, date_array, time_array, year_array, depth_array, lat_array, lng_array = \
            read_regional_catalog(min_mag)

    #   First read the timeseries.txt file.  First determine number of time series
    
    input_file = open('timeseries_EMA_' + LS +'.txt','r')    #   'LS' is either = 'L' or 'S'
    number_ts=0
    for line in input_file:
        number_ts += 1
    input_file.close()
    
    time_bins   = []
    timeseries_EMA = [[] for i in range(number_ts-1)]
    
    input_file = open('timeseries_EMA_' + LS +'.txt','r')
    i=0
    for line in input_file:
        items = line.strip().split()

        if i == 0:
            for j in range(len(items)):
                time_bins.append(float(items[j]))
        else:
            for j in range(len(items)):
                timeseries_EMA[i-1].append(float(items[j]))
                
        i += 1
                
    input_file.close()        

    return time_bins, timeseries_EMA
    
    ######################################################################
    
def get_timeseries_reduced_data():

    mag_array, date_array, time_array, year_array, depth_array, lat_array, lng_array = \
            read_regional_catalog(min_mag)

    #   First read the timeseries.txt file.  First determine number of time series
    
    input_file = open('timeseries_reduced.txt','r')
    number_ts=0
    for line in input_file:
        number_ts += 1
    input_file.close()
    
    time_bins   = []
    timeseries  = [[] for i in range(number_ts-1)]
    
    input_file = open('timeseries_reduced.txt','r')
    i=0
    for line in input_file:
        items = line.strip().split()
#         print('items: ', items)
        if i == 0:
            for j in range(len(items)):
                time_bins.append(float(items[j]))
        else:
            for j in range(len(items)):
                timeseries[i-1].append(float(items[j]))
                
        i += 1
                
    input_file.close()        

    return time_bins, timeseries
    
    ######################################################################
    
def get_yearseries_data():

    input_file = open('yearseries_unbinned.txt','r')
    number_ts=0
    for line in input_file:
        number_ts += 1
    input_file.close()
    
    year_data   =   []
    
    input_file = open('yearseries_unbinned.txt','r')

    for line in input_file:
        items = line.strip().split()
        working_file = []
        for j in range(len(items)):
            working_file.append(float(items[j]))
        year_data.append(working_file)
        
    input_file.close()  
    
    return year_data
    
    ######################################################################
    
def rewrite_catalog(reference_lat,reference_lng):

    input_file = open('USGS_Master.catalog',"r")
    
    augmented_file = open('Augmented_USGS_Catalog.catalog', "w")
    
    for line in input_file:
        items = line.strip().split()
        
        lat1 = float(reference_lat)
        lat2 = float(items[4])
        
        lng1 = float(reference_lng)
        lng2 = float(items[3])
        
        distance_from_ref_point = NTWUtilities.great_circle_distance(lat1, lng1, lat2, lng2)

        items.append(distance_from_ref_point)
        
        print(items[0],items[1],items[2],items[3],items[4],items[5],items[6],items[7], file=augmented_file)
        
    augmented_file.close()
    
    input_file.close()

    return
    
    #################################################################
    
