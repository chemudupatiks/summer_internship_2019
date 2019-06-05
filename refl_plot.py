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
    
    # Plotting the reflectivity image.
    fig = plt.figure()
    plt.imshow(reflectivity.transpose(), cmap = cmap, origin = 'lower')
    plt.colorbar()
    #plt.xticks(times, get_time_hms(times))
    #plt.xlabel('time')
    #plt.ylabel(i)
    #plt.legend()


    # Establishing connection between the figure and the onClick event. 
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()

    # Closing connection. 
    fig.canvas.mpl_disconnect(cid)

def plotXY(x, y):
    fig = plt.figure()
    plt.plot(x, y) 

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show() 

    fig.canvas.mpl_disconnect(cid) 

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
    #max_ref = np.max(reflectivity)
    #min_ref = np.min(reflectivity)

    #plotImage(reflectivity)

    # plot reflectivities of closest upper and lower gate
    altrange = np.array(nc.variables['altrange'][:]) 
    print(altrange.shape) 
    print(alt.shape) 
    print(reflectivity.shape)

    upper_gate_indices = []
    lower_gate_indices = [] 

    for i in range(len(alt)): 
       # print(alt[i])
       temp =  (altrange - alt[i]) 
       temp[temp> 0] = np.nan
       #print(np.where(temp == np.nanmax(temp))[0][0])
       upper_gate.append(reflectivity[i, np.where(temp == np.nanmin(temp))[0][0]])
    #print(np.array(upper_gate).data)
    #print(np.array(upper_gate)[np.where(upper_gate != np.nan)])
    a = reflectivity[:10,:]
    print(a.shape) 
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            print(a[i][j]),
            print( " "), 
        print("\n") 
        
    plotXY(times, upper_gate)
    
    
 

    #print(coordinates) 
    #print(time_alt)

    #crop = reflectivity[int(coordinates[0][0]):int(coordinates[1][0])+1][:]
    #print(crop.shape) 
    #plt.imshow(crop.transpose(), cmap = cmap, origin = 'lower')
    #plt.show()
