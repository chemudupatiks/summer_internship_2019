# author: Krishna Sai Chemudupati
# date  : May 31, 2019 
# email : kchemudu@uwyo.edu

import numpy as np 
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm 
import netCDF4
from sys import argv
import seaborn as sb
import time 
from matplotlib import interactive 
from color_loader import getColorMap as gc

# Global variables 
coordinates = []
time_alt = []
flight_path = [] 
range_cut = 5 # number of nearest gates to consider while searching for the nearest gate with a valid reflectivity value
min_refl = -15
mask_fill_value = -32767

# Figures 
fig = plt.figure()  # refl_image
fig2 = plt.figure() # time series of refl at nearest gates 

# Converts time in seconds from 1970 to time in HH MM SS format (string object)
def get_time_hms(seconds):
    times = []
    for i in seconds:
        times.append(time.strftime('%H:%M:%S', time.gmtime(i)))
    return times

# Click event which recognizes the scroll click to register the x, y values.
def onclick(event):
    if event.button == 2:
       x, y = event.xdata, event.ydata
       try:
          #print(get_time_hms([times[int(x)]]), alt[int(y)])
          print(get_time_hms([times[int(x)+ 3000]]), altrange[int(y)+ 3000])
       except Exception as e:
          #print(e)
          return
       coordinates.append((x, y))
       time_alt.append((get_time_hms([times[int(x)]]), y))


### Plots the reflectivity image
  # Parameters: 
  # reflectivity: Reflectivity matrix
  # threshold   : Minimum value of the reflectivity
  # h_cut       : Cuts the image length (horizontal_cut)
  # v_cut       : Cuts the image height (used it to hide the reflectivities at the terrain-level to better see the colors in the area of interest)
  # x_step      : Step size between each x_tick on the x_axis. Reduce to see more times on the x_axis. 
  # y_step      : Step size between each y_tick on the y_axis. 
def plotImage(reflectivity, threshold, h_cut, v_cut, x_step, y_step): # 3000, 80, 120, 50
    ### Color map configuration
    #cmap = cm.get_cmap('gist_ncar_r') #gist_ncar_r, gist_ncar , and  prism are other colormaps 
    #cmap = gc("../colortables/myBkBlAqGrYeOrReViWh200.gp")
    #cmap = gc("../colortables/BlGrYeOrReVi200.gp")
    cmap = gc("../colortables/BlAqGrYeOrReVi200.rgb")


    ### masking the reflectivity image based on the given threshold. 
    masked_refl = reflectivity # This is causing the rest of the plots to also have masked values according to the threshold.
    masked_refl[np.where(masked_refl < threshold)] = np.nan
    
    ### Plotting 
    fig = plt.figure(1)
    plt.imshow(masked_refl[h_cut:, v_cut:].transpose(), cmap = cmap, origin = 'lower')
    plt.colorbar()
    plt.plot(np.array(flight_path[h_cut:])-v_cut) 
    plt.xticks(np.arange(0, len(times[h_cut:]), x_step), get_time_hms(times[np.arange(h_cut, len(times), x_step)]))
    plt.yticks(np.arange(0, len(altrange[v_cut:]), y_step), altrange[np.arange(v_cut, len(altrange), y_step)])
    plt.xlabel("Time (UTC)") 
    plt.ylabel("Altitude (m MSL)")

    
    # Establishing connection between the figure and the onClick event. 
    cid = fig.canvas.mpl_connect('button_press_event', onclick)


### Plots the reflectivity vs time at the nearest gates with some reflectivity
  # Parameters: 
  # r_c : Number of gates above and below the flight level being searched for existing refl value
  # threshold : The min refl value, all values below it are masked
  # x_step    : Step size between two x_ticks, decrease value to increase the number of x_ticks on the x axis
def plotNearestGateReflectivity(r_c, threshold, x_step):
    upper_gate = []
    lower_gate = [] 

    for i in range(len(alt)):
       abs_diff = np.abs(altrange - alt[i]) 
       flight_level = np.argmin(abs_diff) 
       refl_all_gates = reflectivity.data[i,:]
       refl_all_gates[:(flight_level-range_cut)] = np.nan 
       refl_all_gates[(flight_level+range_cut):] = np.nan
       refl_all_gates[refl_all_gates == mask_fill_value] = np.nan
       non_nan_refl_index = np.array(np.where(~(np.isnan(refl_all_gates)))).flatten()
       diff_gates = non_nan_refl_index - flight_level
       try: 
           upper_gate.append(reflectivity[i][non_nan_refl_index[np.where(diff_gates == (np.min(diff_gates[diff_gates > 0])))[0][0]]])
       except Exception as e:
           #print("upper_gate: ", e) 
           upper_gate.append(np.nan) 
       try: 
           lower_gate.append(reflectivity[i][non_nan_refl_index[np.where(diff_gates == (np.max(diff_gates[diff_gates < 0])))[0][0]]])
       except Exception as e:
           #print("lower_gate: ", e) 
           lower_gate.append(np.nan) 
    
    #print(len(lower_gate)) 
    #print(len(upper_gate)) 
    #plotXY(times, upper_gate)
    #plotXY(times, lower_gate)
    
    upper_gate = np.array(upper_gate) 
    lower_gate = np.array(lower_gate)   
    if(threshold is not None):  
      upper_gate[np.where(np.isnan(upper_gate))] = mask_fill_value 
      lower_gate[np.where(np.isnan(lower_gate))] = mask_fill_value 
      lower_gate[lower_gate < threshold] = np.nan
      upper_gate[upper_gate < threshold] = np.nan
     
    #print(np.where(~(np.isnan(upper_gate))))
     
    fig2 = plt.figure(2)
    plt.plot(times, upper_gate, label="Closest Gate above flight track") 
    plt.plot(times, lower_gate, label="Closest Gate below flight track")
    plt.legend()
    plt.xlabel("Time (UTC)") 
    plt.ylabel("Reflectivity (dBZ)")
    plt.xticks(times[np.arange(0, len(times), x_step)], get_time_hms(times[np.arange(0, len(times), x_step)]))
    cid = fig2.canvas.mpl_connect('button_press_event', onclick)

### Main script ###
n_args = len(argv)
if n_args <= 1:
    print("please provide the filename(s) to be parsed") 
    quit()
else:
    filename = argv[1]
    try:
        # Opens the netCdf file and exits if incorrect filename  is provided. 
        nc = netCDF4.Dataset(filename, 'r')
    except:
        print("Please check the validity of the input files")
        quit()

    # extracts reflectivities, times, and altitude from the cdf file. 
    reflectivity  = nc.variables['reflectivity'][:]
    times = nc.variables['time'][:]
    alt = nc.variables['ALT'][:]
    altrange = np.array(nc.variables['altrange'][:]) 

    ### Filling the flight_path array (global)
    for i in range(len(alt)):
       abs_diff = np.abs(altrange - alt[i]) 
       flight_path.append(np.argmin(abs_diff))

    # plot refl image
    plotImage(reflectivity, min_refl, 3000, 80, 120, 50)

    # plot reflectivities of closest upper and lower gate
    plotNearestGateReflectivity(range_cut, min_refl, 50)  

    # Show all plots   
    plt.show() 

    # Close all connections
    fig2.canvas.mpl_disconnect(cid)
    fig.canvas.mpl_disconnect(cid) 
