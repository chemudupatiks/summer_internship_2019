# author: Krishna Sai Chemudupati
# date  : July 2, 2019 
# email : kchemudu@uwyo.edu

import numpy as np 
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm 
import netCDF4
from sys import argv
import time
from matplotlib import interactive 
from color_loader import getColorMap as gc
from datetime import timedelta, datetime 



def XY_scatter_plot(filename, x, x_time, y, y_time, date = None, s_time = None, e_time = None):    
    x_st_index = -1; x_end_index = -1; y_st_index = -1; y_end_index = -1
    if date and s_time and e_time:
        dt = datetime(1970, 1, 1)  
        start_time = date +" "+ s_time 
        end_time = date +" "+ e_time 
        start_time = (datetime.fromtimestamp(time.mktime(time.strptime(start_time, "%m/%d/%Y %H.%M.%S")))-dt).total_seconds()
        end_time = (datetime.fromtimestamp(time.mktime(time.strptime(end_time, "%m/%d/%Y %H.%M.%S")))-dt).total_seconds()
        x_st_index = np.where(x_time.astype(int) == start_time)[0][0]
        x_end_index = np.where(x_time.astype(int) == end_time)[0][-1]
        y_st_index = np.where(y_time.astype(int) == start_time)[0][0]
        y_end_index = np.where(y_time.astype(int) == end_time)[0][-1]

    gcs_csv = pd.read_csv("../csv/gcs.csv")
    for i in range(gcs_csv.shape[0]):
        print(gcs_csv.iloc[0]) 



    
def get_time_intervals(t, threshold, h_cut, upper):
    if upper:
        time_in_cells = np.where(upper_gate >= threshold)[0]
    else:
        time_in_cells = np.where(lower_gate >= threshold)[0]

    time_in_cells = time_in_cells[time_in_cells >= h_cut]
    time_in_cells = t[time_in_cells]

    return time_in_cells


n_args = len(argv)
if n_args <= 1:
    print("please provide the filename(s) to be parsed") 
    quit()
else:
    filename = argv[1]
    filename2 = argv[2]
    try:
        # Opens the netCdf file and exits if incorrect filename  is provided. 
        nc = netCDF4.Dataset(filename, 'r')
        nc_sensor = netCDF4.Dataset(filename2, 'r') 
    except:
        print("Please check the validity of the input files")
        quit()

    # extracts reflectivities, times, and altitude from the cdf file. 
    reflectivity  = nc.variables['reflectivity'][:]
    times = nc.variables['time'][:]
    altitude = nc.variables['ALT'][:]
    alt_range = np.array(nc.variables['altrange'][:])
    vel_corr = nc.variables['velocity_corrected'][:]
    
    full_time = nc_sensor.variables['time'][:]
    dt = datetime(2017, 1, 1) 
    dt2 = datetime(1970, 1, 1) 
    full_time += (dt - dt2).total_seconds() 
    full_time = np.array(full_time).astype(int)
    rosemount = nc_sensor.variables['rlwc']
    dmt100 = nc_sensor.variables['lwc100']
    dmtcdp = nc_sensor.variables['cdplwc_NRB']
    nevz = nc_sensor.variables['nevlwc']
    total_nevz = nc_sensor.variables['nevtwc']
    eddy = nc_sensor.variables['turb']
    flight_path  = []

    #print(vel_corr.shape) 

    ### Filling the flight_path array (global)
    for i in range(len(altitude)):
       abs_diff = np.abs(alt_range - altitude[i]) 
       flight_path.append(np.argmin(abs_diff))


    upper_gate = []
    lower_gate = []
    mask_fill_value = -32767
    r_c = 5
    threshold = -8 

    for i in range(len(altitude)):
       abs_diff = np.abs(alt_range - altitude[i]) 
       flight_level = np.argmin(abs_diff) 
       refl_all_gates = reflectivity.data[i,:]
       refl_all_gates[:(flight_level-r_c)] = np.nan 
       refl_all_gates[(flight_level+r_c):] = np.nan
       refl_all_gates[refl_all_gates == mask_fill_value] = np.nan
       non_nan_refl_index = np.array(np.where(~(np.isnan(refl_all_gates)))).flatten()
       diff_gates = non_nan_refl_index - flight_level
       try: 
           upper_gate.append(refl[i][non_nan_refl_index[np.where(diff_gates == (np.min(diff_gates[diff_gates > 0])))[0][0]]])
       except Exception as e:
           #print("upper_gate: ", e) 
           upper_gate.append(np.nan) 
       try: 
           lower_gate.append(refl[i][non_nan_refl_index[np.where(diff_gates == (np.max(diff_gates[diff_gates < 0])))[0][0]]])
       except Exception as e:
           #print("lower_gate: ", e) 
           lower_gate.append(np.nan) 
    
    
    upper_gate = np.array(upper_gate) 
    lower_gate = np.array(lower_gate)


    if(threshold is not None):  
      upper_gate[np.where(np.isnan(upper_gate))] = mask_fill_value 
      lower_gate[np.where(np.isnan(lower_gate))] = mask_fill_value 
    
      lower_gate[lower_gate < threshold] = mask_fill_value
      upper_gate[upper_gate < threshold] = mask_fill_value 


    XY_scatter_plot(filename, nevz, full_time, vel_corr, times) 



