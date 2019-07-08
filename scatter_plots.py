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

def expand_array(x_sub, x_sub_time, y_sub, y_sub_time):
    low_res = []; low_res_time = []; high_res= []; high_res_time = []; low_res_exp = []; low_res_time_exp = []
    if x_sub.shape[0] > y_sub.shape[0]:
        low_res = y_sub; low_res_time = y_sub_time 
        high_res = x_sub; high_res_time = x_sub_time 
    else: 
        low_res = x_sub; low_res_time = x_sub_time 
        high_res = y_sub; high_res_time = y_sub_time
    for i in range(low_res.shape[0]): 
        time_match = np.where(high_res_time.astype(int) == low_res_time[i]) 
        #print(time_match[0]) 
        for j in range(time_match[0].shape[0]): 
            low_res_time_exp.append(low_res_time[i]) 
            low_res_exp.append(low_res[i]) 
    #print(np.array(low_res_exp).shape) 
    #print(high_res.shape)
    #print(low_res_time_exp[-5:]) 
    #print(high_res_time[-5:])
    return np.array(low_res_exp), np.array(low_res_time_exp) 



# Have to pass time arrays with seconds from Jan 1 1970. Also the x, and y arrays must be linear (1D).
def XY_scatter_plot(filename, x, x_time, y, y_time, date = None, s_time = None, e_time = None, xlabel = None, ylabel = None, title = None):    
    x_st_index = -1; x_end_index = -1; y_st_index = -1; y_end_index = -1
    x_sub_time = x_time[:]
    y_sub_time = y_time[:]
    x_sub = x[:] 
    y_sub = y[:] 
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
        x_sub_time = x_time[x_st_index:x_end_index+1]
        y_sub_time = y_time[y_st_index:y_end_index+1]
        x_sub = x[x_st_index:x_end_index+1]
        y_sub = y[y_st_index:y_end_index+1] 
        #print(x_sub_time.shape)
        #print(x_sub.shape) 
        #print(y_sub.shape) 
        #print(y_sub_time.shape)
        #print(start_time) 
        #print(end_time) 
        #print(x_sub_time[-5:]) 
        #print(y_sub_time[-5:]) 
    if(x_sub_time.shape[0] >  y_sub_time.shape[0]): 
        y_sub, y_sub_time = expand_array(x_sub, x_sub_time, y_sub, y_sub_time) 
    elif y_sub_time.shape[0] > x_sub_time.shape[0]:
        x_sub, x_sub_time = expand_array(x_sub, x_sub_time, y_sub, y_sub_time) 
    #print(x_sub_time.shape)
    #print(x_sub.shape) 
    #print(y_sub.shape) 
    #print(y_sub_time.shape)
    inside_gcs = np.array([False] * x_sub_time.shape[0]) 
    #print(color_scheme)

    gcs_csv = pd.read_csv("../csv/gcs.csv")
    for i in range(gcs_csv.shape[0]):
        #print(gcs_csv.iloc[i])
        #print("\n\n\n") 
        array = np.arange(int(gcs_csv.loc[i]['start_seconds']) ,int(gcs_csv.loc[i]['end_seconds'])+1)
        for e in array:
            if e >= start_time and e <= end_time: 
                x_index = np.where(x_sub_time == e)[0]
                #print(x_index)
                inside_gcs[x_index] = True
                #y_index = np.where(y_sub_time == e) 
                #print(x_index) 
                #print(y_index) 
    #print(color_scheme)
    #plt.scatter(x_sub, y_sub, c = color_scheme)
    x_sub[x_sub == mask_fill_value] = np.nan
    y_sub[y_sub == mask_fill_value] = np.nan 

    plt.plot(
            x_sub[inside_gcs], y_sub[inside_gcs], 'bo', 
            x_sub[~inside_gcs], y_sub[~inside_gcs], 'ro', 
    )
    if xlabel is not None: plt.xlabel(xlabel) 
    if ylabel is not None: plt.ylabel(ylabel) 
    if title is not None : plt.title(title)
    plt.legend(['inside gcs','outside gcs'])
    plt.grid(True)
    plt.show()

    
def get_time_intervals(t, threshold, h_cut, upper_gate, lower_gate, upper):
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
    
    full_time = np.array(nc_sensor.variables['time'][:]) 
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


    upper_gate_refl = []
    lower_gate_refl = []
    upper_gate_vel = []
    lower_gate_vel= []
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
           upper_gate_refl.append(reflectivity[i][non_nan_refl_index[np.where(diff_gates == (np.min(diff_gates[diff_gates > 0])))[0][0]]])
           upper_gate_vel.append(vel_corr[i][non_nan_refl_index[np.where(diff_gates == (np.min(diff_gates[diff_gates > 0])))[0][0]]])
       except Exception as e:
           #print("upper_gate: ", e) 
           upper_gate_refl.append(np.nan)
           upper_gate_vel.append(np.nan) 
       try: 
           lower_gate_refl.append(reflectivity[i][non_nan_refl_index[np.where(diff_gates == (np.max(diff_gates[diff_gates < 0])))[0][0]]])
           lower_gate_vel.append(vel_corr[i][non_nan_refl_index[np.where(diff_gates == (np.max(diff_gates[diff_gates < 0])))[0][0]]])
       except Exception as e:
           #print("lower_gate: ", e) 
           lower_gate_refl.append(np.nan) 
           lower_gate_vel.append(np.nan) 
    
    
    upper_gate_refl = np.array(upper_gate_refl) 
    lower_gate_refl = np.array(lower_gate_refl)
    upper_gate_vel = np.array(upper_gate_vel) 
    lower_gate_vel = np.array(lower_gate_vel)


    if(threshold is not None):  
      upper_gate_refl[np.where(np.isnan(upper_gate_refl))] = mask_fill_value 
      lower_gate_refl[np.where(np.isnan(lower_gate_refl))] = mask_fill_value
   
      lower_gate_refl[lower_gate_refl < threshold] = mask_fill_value
      upper_gate_refl[upper_gate_refl < threshold] = mask_fill_value 
      lower_gate_vel[lower_gate_refl < threshold] = mask_fill_value
      upper_gate_vel[upper_gate_refl < threshold] = mask_fill_value 

    #print(lower_gate_vel[lower_gate_refl>= threshold]) 
    #print(nevz) 


    #XY_scatter_plot(filename, nevz, full_time, vel_corr, times.astype(int))
    #XY_scatter_plot(filename, nevz, full_time, lower_gate_vel, times.astype(int), date="03/09/2017", s_time = "14.23.11", e_time="14.24.38", xlabel = "Nevzorov LWC", ylabel = "lower gate Doppler vel") 
    #XY_scatter_plot(filename, lower_gate_vel, times.astype(int),  nevz, full_time, date="03/09/2017", s_time = "14.23.11", e_time="14.24.38", xlabel = "Nevzorov LWC", ylabel = "lower gate Doppler vel") 
    XY_scatter_plot(filename, nevz, full_time, upper_gate_vel, times.astype(int), date="03/09/2017", s_time = "14.23.11", e_time="14.24.38", xlabel = "Nevzorov LWC", ylabel = "upper gate Doppler vel")

