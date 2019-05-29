import numpy as np 
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm 
import netCDF4
from sys import argv
import seaborn as sb

def get_time_hms(seconds):
    pass

n_args = len(argv)
if n_args <= 1:
    print("please provide the filename(s) to be parsed") 
    quit()
else:
    filename = argv[1]
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
    time = nc.variables['time'][:]
    max_ref = np.max(reflectivity)
    min_ref = np.min(reflectivity)
    cmap = cm.get_cmap('inferno')
    ext_time  = []
    for i in range(len(time)):
        ext_time.append(np.ones(len(reflectivity[0]), dtype = 'float')*time[i])
    ext_time  = np.array(ext_time).flatten()
    print(len(reflectivity.flatten()))
    print(len(np.array(ext_time).flatten()))
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
    plt.plot(time, nc.variables['ALT'][:])
    plt.show()
