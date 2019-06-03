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

   
    print(coordinates) 
    print(time_alt)
