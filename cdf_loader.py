import numpy as np 
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm 
import netCDF4
from sys import argv
import seaborn as sb
import time 

def get_time_hms(seconds):
    times = []
    for i in seconds:
        times.append(time.strftime('%H:%M:%S', time.gmtime(i%86399)))
    return times

n_args = len(argv)
if n_args <= 1:
    print("please provide the filename(s) to be parsed") 
    quit()
else:
    filename = argv[1]
    var_plot = argv[2:]
    try:
        nc = netCDF4.Dataset(filename, 'r')
    except:
        print("Please check the validity of the input files")
        quit()
    #print(nc.groups.values)
    #print(nc.dimensions)
    #print(nc.variables)
    #print(nc.variables.keys())
    #print(nc.variables['reflectivity'])
    #print(nc.variables['ALT'])
    reflectivity  = nc.variables['reflectivity'][:]
    times = nc.variables['time'][:]
    #max_ref = np.max(reflectivity)
    #min_ref = np.min(reflectivity)
    cmap = cm.get_cmap('gist_ncar_r') #gist_ncar_r, prism, 
    #ext_time  = []
    #for i in range(len(time)):
    #    ext_time.append(np.ones(len(reflectivity[0]), dtype = 'float')*time[i])
    #ext_time  = np.array(ext_time).flatten()
    #print(len(reflectivity.flatten()))
    #print(len(np.array(ext_time).flatten()))
    #print(len(cmap.colors))
    #print(nc.variables['time'][:])
    #plt.plot(nc.variables['time'][:], nc.variables['nevlwc'][:])
    #plt.plot(nc.variables['time'][:], nc.variables['nevtwc'][:])
    #plt.show()
    #plt.plot(time, reflectivity)
    #plt.show()
    #norm = matplotlib.colors.Normalize(vmin = min_ref, vmax = max_ref)
    #plt.scatter(ext_time, reflectivity.flatten(), c = reflectivity.flatten(), vmin = min_ref, vmax = max_ref)
    #heat_map = sb.heatmap(np.array(reflectivity), cmap = cmap)
    #plt.plot(nc.variables['time'][:], nc.variables['nevtwc'][:])
    #plt.plot(time, nc.variables['ALT'][:])
    #plt.show()
    #print(get_time_hms(times)[-10:]) 
    for i in var_plot:
        #plt.plot(times, nc.variables[i][:])
        plt.imshow(reflectivity.transpose(), cmap = cmap, origin = 'lower')
        plt.colorbar()
        #plt.xticks(times, get_time_hms(times))
        #plt.xlabel('time')
        #plt.ylabel(i)
        #plt.legend()
        plt.show()
    quit() 
