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

# Global variables to hold the selected times and y-value( gate?) 
coordinates = []
time_alt = []

# Converts time in seconds from 1970 to time in HH MM SS format  
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
          print(get_time_hms([times[int(x)]]), y)
       except Exception as e:
          #print(e)
          return
       coordinates.append((x, y))
       time_alt.append((get_time_hms([times[int(x)]]), y))

def plotImage(reflectivity):
    # Color map configuration
    cmap = cm.get_cmap('gist_ncar_r') #gist_ncar_r, gist_ncar , and  prism are other colormaps 
    #cmap = gc("../colortables/myBkBlAqGrYeOrReViWh200.gp")
    #cmap = gc("../colortables/BlGrYeOrReVi200.gp")
    cmap = gc("../colortables/BlAqGrYeOrReVi200.rgb")
    # Plotting the reflectivity image
    masked_refl = reflectivity
    #masked_refl = nc.variables['reflectivity'][:] 
    #masked_refl[np.where(masked_refl < -15)] = np.nan
    #print(np.where(masked_refl < -15))
    #quit()
    fig = plt.figure(1)
    plt.imshow(masked_refl[3000:, 80:].transpose(), cmap = cmap, origin = 'lower')
    plt.colorbar()
    #plt.xticks(times, get_time_hms(times))
    #print(np.arange(0, len(times), 120), get_time_hms(times[np.arange(0, len(times), 120)]))
    plt.xticks(np.arange(0, len(times[3000:]), 120), get_time_hms(times[np.arange(3000, len(times), 120)]))
    plt.yticks(np.arange(0, len(altrange[80:]), 50), altrange[np.arange(80, len(altrange), 50)])
    #plt.xlabel('time')
    #plt.ylabel(i)
    #plt.legend()

    
    # Establishing connection between the figure and the onClick event. 
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    #interactive(True)
    plt.show()

    # Closing connection. 
    fig.canvas.mpl_disconnect(cid)

def plotXY(x, y):
    fig = plt.figure(2)
    plt.plot(x, y) 

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    #plt.show() 

    fig.canvas.mpl_disconnect(cid) 

### Main script ###
interactive(True)
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

    #max_ref = np.max(reflectivity)
    #min_ref = np.min(reflectivity)

    plotImage(reflectivity)

    # plot reflectivities of closest upper and lower gate
    altrange = np.array(nc.variables['altrange'][:]) 
    #print(altrange.shape) 
    #print(alt.shape) 
    #print(reflectivity.shape)

    upper_gate = []
    lower_gate = [] 

    range_cut = 5 
    min_refl = -10
    mask_fill_value = -32767
    for i in range(len(alt)):
       abs_diff = np.abs(altrange - alt[i]) 
       flight_level = np.argmin(abs_diff) 
       refl_all_gates = reflectivity.data[i,:]
       refl_all_gates[:(flight_level-range_cut)] = np.nan 
       refl_all_gates[(flight_level+range_cut):] = np.nan
       refl_all_gates[refl_all_gates == mask_fill_value] = np.nan
       non_nan_refl_index = np.array(np.where(~(np.isnan(refl_all_gates)))).flatten()
       #print(non_nan_refl_index)
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

       #print(upper_gate_indices) 
       #print(lower_gate_indices)
       #print(non_nan_refl_index)  
       #print(refl_all_gates) 
       #print(flight_level) 
       #print(abs_diff) 
       #quit()
    
    #print(len(lower_gate)) 
    #print(len(upper_gate)) 
    #quit()

    #plotXY(times, upper_gate)
    #plotXY(times, lower_gate)
    temp  = lower_gate[:] 
    temp2 = upper_gate[:]
    
    upper_gate = np.array(upper_gate) 
    lower_gate = np.array(lower_gate) 
    
    #upper_gate[np.where(np.isnan(upper_gate))] = mask_fill_value 
    #lower_gate[np.where(np.isnan(lower_gate))] = mask_fill_value 

    #lower_gate[lower_gate < min_refl] = np.nan
    #upper_gate[upper_gate < min_refl] = np.nan
     
    #print(np.where(~(np.isnan(upper_gate))))
     
    plt.figure(2)
    plt.plot(times, upper_gate, label="Closest Gate above flight track") 
    plt.plot(times, lower_gate, label="Closest Gate below flight track")
    plt.legend()
    plt.xticks(times[np.arange(0, len(times), 120)], get_time_hms(times[np.arange(0, len(times), 120)]))
    interactive(False)
    plt.show()
