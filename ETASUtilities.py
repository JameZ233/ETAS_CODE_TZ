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
import textwrap
import matplotlib.mlab as mlab
#from matplotlib.pyplot import figure, show
from matplotlib import cbook
import numpy as np
from numpy import *
# from mpl_toolkits.basemap import Basemap
from array import array
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.patches as patches

import datetime
import dateutil.parser
from datetime import timedelta, datetime

import urllib.request, urllib.parse, urllib.error
import urllib.request, urllib.error, urllib.parse
import os

import math
from math import exp
from math import log

from matplotlib import cm

import http.client
from urllib.error import HTTPError

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

import random
import operator

import requests
from itertools import product

from osgeo import ogr, osr

    ######################################################################
    
def convert_partial_year(number):

    year = int(number)
    d = timedelta(days=(number - year)*(365))
    day_one = datetime(year,1,1)
    date = d + day_one
    
    form_date_dash = str(date.date())
    form_date = form_date_dash.split('-')
    form_date = str(form_date[1]) + '/' + str(form_date[2]) + '/' + str(form_date[0])
    
    return form_date_dash, form_date
    
    #   .................................................................

def mean_val(number_array):

    #   The obvious: returns the mean of an array of numbers

    N_Dim_mean = len(number_array)
    mean = sum(number_array)/float(N_Dim_mean)

    return mean

    #   .................................................................

def std_var_val(number_array):

    #   The obvious: returns the standard deviation of an array of numbers

    N_Dim_var = len(number_array)
    mean_of_array = mean_val(number_array)

    adjusted_arguments = [((i-mean_of_array)**2)/float(N_Dim_var-1) for i in number_array]   # Sample variance

    variance = sum(adjusted_arguments)   

    standard_deviation = math.sqrt(variance) 

    return (standard_deviation, variance)
    
    #   .................................................................
    
def median_value_intervals(number_array,count,median_flag):

    future_array   =   []
    
    total_number = len(number_array)    #   Number of intervals
    
    if median_flag == 'true':                #   Use only positive future intervals
        for index in range(0,total_number):
            current_number = number_array[index]
            if current_number >= count:
                future_count = current_number - count
                future_array.append(future_count)
                
    if median_flag == 'false':                #   Use all future intervals
        for index in range(0,total_number):
            current_number = number_array[index]
            future_count = current_number - count
            future_array.append(future_count)

        number_usable_intervals = len(future_array)           

    number_usable_intervals = len(future_array)
    
    mean_count = '(No Data)'
    mean_stdev = '(No Data)'
    median_count = '(No Data)'
    percent_25 = '(No Data)'
    percent_75 = '(No Data)'
    percent_99 = '(No Data)'
        
    if number_usable_intervals > 2:
        mean_count = str(int(round(float(mean_val(future_array)),0)))
        mean_stdev, mean_var = std_var_val(future_array)
        median_count = str(int(np.median(future_array)))
        
        #  In the next statement, notice we need to be aware of what 25%, 75% and 99% really means 
        #       Or in other words, since are dealing with survivor or exceedance probabilities,
        #       it means that only a few (e.g., 25%) of intervals are long, and most of them (e.g. 75%) are short
        #       etc.
        
        percent_25, percent_75, percent_99= np.percentile(future_array,q=[75.0,25.0,1.0])
        percent_25 = str(int(percent_25))
        percent_75 = str(int(percent_75))
        percent_99 = str(int(percent_99))
        
    return (median_count,percent_25,percent_75,percent_99,mean_count,mean_stdev,number_usable_intervals)

    #   .................................................................

def linfit(x, y, n):

#
#     This program fits a straight line to N data.  The data are
#     assumed to be error-free.  
#
#     Definitions:
#
#     x[i]:  Abscissas of the data
#     y[i]:  Ordinates of the data
#     slope: The resulting best fit slope
#     cept:  The resulting best fit y-intercept
#     errs:  The standard error for the slope
#     errc:  The standard error for the intercept
#     n:     The exact number of data to fit
#
#
#   NOTE!!:     The *exact number* n of data to fit *must equal* = len(x) = len(y)
#               Otherwise, the code will blow up
#
#
#    n = int(len(x))

    n = int(n)

    ata = [zeros(2),zeros(2)]
    aty = [zeros(2),zeros(2)]
    atainv = [zeros(2),zeros(2)]
#    
    sumx2 = 0.
    xsum = 0.
    yxsum = 0.
    ysum = 0.

    for i in range(0,n):
        sumx2 = sumx2 + x[i] * x[i]
        xsum = xsum + x[i]
        yxsum = yxsum + y[i] * x[i]
        ysum = ysum + y[i]
#
#    ata[0][0] = sumx2
#    ata[0][1] = xsum
#    ata[1][0] = xsum
#    ata[1][1] = float(n)

    ata[0][0] = sumx2
    ata[1][0] = xsum
    ata[0][1] = xsum
    ata[1][1] = float(n)
#
    aty[0] = yxsum
    aty[1] = ysum

#
    det = ata[0][0] * ata[1][1] - ata[0][1] * ata[1][0]
    atainv[0][0] = ata[1][1]/det
    atainv[0][1] = -ata[0][1]/det
    atainv[1][0] = -ata[1][0]/det
    atainv[1][1] = ata[0][0]/det
#
    slope = atainv[0][0] * aty[0] + atainv[0][1] * aty[1]
    cept = atainv[1][0] * aty[0] + atainv[1][1] * aty[1]

    s2 = 0
    for i in range(0,n):
        s2 = s2 + (y[i] - cept - slope * x[i])**2

    s2 = s2 / (float(n) - 2.)
      
    errs = math.sqrt( float(n) * s2 / det )
    errc = math.sqrt( s2 * sumx2 / det)

#   print slope, cept, errs, errc, s2

    return (slope, cept, errs, errc, s2)

    #   Usage in calling program:  slope, cept, errs, errc, s2 = linfit(x,y)

    ##################################################################


def change_in_latitude(km):

    #   http://www.johndcook.com/blog/2009/04/27/converting-miles-to-degrees-longitude-or-latitude/

    # Distances are measured in km.
    # Longitudes and latitudes are measured in degrees.
    # Earth is assumed to be perfectly spherical.

    earth_radius = 6371.0
    degrees_to_radians = math.pi/180.0
    radians_to_degrees = 180.0/math.pi

    #Given a distance north, return the change in latitude in degrees.

    return (km/earth_radius)*radians_to_degrees

    ##################################################################

def change_in_longitude(latitude, km):

    #   http://www.johndcook.com/blog/2009/04/27/converting-miles-to-degrees-longitude-or-latitude/

    # Distances are measured in km.
    # Longitudes and latitudes are measured in degrees.
    # Earth is assumed to be perfectly spherical.

    earth_radius = 6371.0
    degrees_to_radians = math.pi/180.0
    radians_to_degrees = 180.0/math.pi

    # Given a latitude and a distance west, return the change in longitude in degrees.
    # Find the radius of a circle around the earth at given latitude.

    r = earth_radius*math.cos(latitude*degrees_to_radians)
    return (km/r)*radians_to_degrees

    ##################################################################

def draw_big_circle(CircleCenterLat, CircleCenterLng, CircleRadius):

    # Build an array of x-y values in degrees defining a circle of the required radius

    twopi = 6.283185306

    number_points   = 100
    delta_angle     = twopi/float(number_points)

    x_circle_dg      = np.zeros(number_points)
    y_circle_dg      = np.zeros(number_points)

    angle = 0.0     #   start at due north

    for i in range(0,number_points):
        x_comp = CircleRadius * math.sin(angle) 
        y_comp = CircleRadius * math.cos(angle) 

        y_circle_dg[i] = CircleCenterLat + change_in_latitude(y_comp)
        x_circle_dg[i] = CircleCenterLng + change_in_longitude(CircleCenterLat, x_comp)

        angle += delta_angle

    return (x_circle_dg, y_circle_dg)

    ##################################################################
    
def compute_great_circle_distance(lat_1, lng_1, lat_2, lng_2):

    # Build an array of x-y values in degrees defining a circle of the required radius

    pic = 3.1415926535/180.0
    Radius = 6371.0

    lat_1 = float(lat_1) * pic
    lng_1 = float(lng_1) * pic
    lat_2 = float(lat_2) * pic
    lng_2 = float(lng_2) * pic
    
    delta_lng = lng_1 - lng_2
    
    delta_radians = math.sin(lat_1)*math.sin(lat_2) + math.cos(lat_1)*math.cos(lat_2)*cos(delta_lng)
    if delta_radians > 1.0:
        delta_radians = 1.0
    delta_radians = math.acos(delta_radians)
    
    great_circle_distance = delta_radians * Radius

    return great_circle_distance

    ##################################################################

def mesh_seismicity(Country, NELat, NELng, SWLat, SWLng, MagLo):

    #   Use 0.1 degree coarse graining in space and 1/100 year in time.
    #   Start by finding number of space-time grid points
    
    delta_lat = float(NELat) - float(SWLat)
    delta_lng = float(NELng) - float(SWLng)
    
    down_scale = 1.0  #   Number of finer scale boxes per degree
    
    n_lat = int(delta_lat*down_scale)
    n_lng = int(delta_lng*down_scale)
    
    n_boxes = n_lat*n_lng
    
    n_times = 2000      #   Arbitrary
    
    #   Now open the working file with the earthquake data

    with open("USGS_master.catalog") as f:
        data_file = f.readlines()

    lat_list    = []
    lng_list    = []
    mag_list    = []
    time_list   = []
    
    filter_magnitude = 5.0
    filter_year = 1980.0

#    if Country == 'USA':
#        filter_magnitude = 4.75
        
    for line in data_file:
        items = line.strip().split()
        if float(items[5]) >= filter_magnitude and float(items[2]) >= filter_year:
            lat_list.append(items[4])
            lng_list.append(items[3])
            mag_list.append(items[5])
            time_list.append(items[2])
        
    n_eqs = len(time_list)
    
    time_first_eq = float(time_list[0])
    time_last_eq  = float(time_list[n_eqs-1])
    
    delta_time = (time_last_eq - time_first_eq)
    
    # We need to build an 2D array of events in lat-lng boxes and coarse grained time
    
    data_array_full = np.zeros((n_boxes,n_times))    #   All spatial boxes are filled with zeros, i.e., no events, initially
    nonzero_boxes = np.zeros(n_boxes)
    
    for i in range(0,n_eqs):

        latitude = float(float(NELat)-float(lat_list[i]))
        longitude= float(float(lng_list[i]) - float(SWLng))
        
        lat_index = int((latitude/delta_lat) * float(n_lat))
        lng_index = int((longitude/delta_lng) * float(n_lng))
        
        box_index = int(n_lng * lat_index + lng_index)       #   box_index goes left to right (lng), top to bottom (lat)
        
        time_index = int(((float(time_list[i]) - time_first_eq) / delta_time)*n_times)        #   Equidistantly tabulated times??
        
        if time_index == n_times:
            time_index = time_index - 1
        
        data_array_full[box_index][time_index] += 1.0
        
        nonzero_boxes[box_index] = 1    #   Only the nonzero boxes will be considered
        
    total_nonzero_boxes = int(sum(nonzero_boxes))
    
    #   Redefine n_boxes
    n_boxes_full = n_boxes
    n_boxes = total_nonzero_boxes
    
    data_array = np.zeros((n_boxes,n_times))
    
    #   Redefine data_array to only contain boxes with events
    i=-1
    for box_index in range(0,n_boxes_full):
        if nonzero_boxes[box_index] == 1:
            i += 1
            for time_index in range(0,n_times):
                data_array[i][time_index] = data_array_full[box_index][time_index]

    return (data_array, n_times, n_boxes, time_first_eq, time_last_eq, down_scale, filter_magnitude, filter_year)
    
    ##################################################################
    
def createCircleAroundWithRadius(lat, lon, radiusKm):

    #   This method draws circles with Haversine functions
    #       Inputs of lat lon in degrees

    ring = ogr.Geometry(ogr.wkbLinearRing)
    
    latArray = []
    lonArray = []
    
    for brng in range(0,366,6):
    
        lat2, lon2 = getLocation(lat,lon,brng,radiusKm)
        latArray.append(lat2)
        lonArray.append(lon2)
        
    return lonArray,latArray

    ##################################################################

def getLocation(lat1, lon1, brng, distanceKm):

    lat1 = float(lat1) * math.pi/ 180.0
    lon1 = float(lon1) * math.pi / 180.0
    #earth radius
    R = 6371.0
    #R = ~ 3959 MilesR = 3959

    distanceKm = distanceKm/R

    brng = (float(brng) / 90.0)* math.pi / 2.0

    lat2 = math.asin(math.sin(lat1) * math.cos(distanceKm) + math.cos(lat1) * math.sin(distanceKm) * math.cos(brng))
    lon2 = lon1 + math.atan2(math.sin(brng)*math.sin(distanceKm)*math.cos(lat1),math.cos(distanceKm)-math.sin(lat1)*math.sin(lat2))

    lon2 = 180.0 * lon2/ math.pi
    lat2 = 180.0 * lat2/ math.pi
    
#    print lat2, lon2

    return lat2, lon2
    
    ##################################################################
    
def calc_information(y_list):

    N = sum(y_list)
    
    if int(N) > 0 and int(N) != len(y_list):       #   Not all zeros or ones
        pr = float(N)/float(len(y_list))
        qr = 1. - pr
        shannon_info = - (  pr * math.log(pr,2) + qr * math.log(qr,2) )
    else:
        shannon_info = 0.0

    return shannon_info
    
    ##################################################################
    