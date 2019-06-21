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
import time
from matplotlib import interactive 
from color_loader import getColorMap as gc
from datetime import timedelta, datetime 

# Global variables 
coordinates = []
time_alt = []
flight_path = [] 
range_cut = 5 # number of nearest gates to consider while searching for the nearest gate with a valid reflectivity value
min_refl = None
min_vel = -5
mask_fill_value = -32767
hor_cut = 3000 
ver_cut = 80
refl_upper_gate = np.array([])
refl_lower_gate = np.array([])
masked_times = np.array([])
# Figures 
fig = plt.figure()  # refl_image
cid = 0
fig2 = plt.figure() # time series of refl at nearest gates 
cid2 = 0
fig3 = plt.figure() # time series of refl at nearest gates 
cid3 = 0


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
          print(get_time_hms([times[int(x)+ hor_cut]]), alt_range[int(y)+ ver_cut])
       except Exception as e:
          print(e)
          return
       coordinates.append((x + hor_cut, y + ver_cut))
       time_alt.append((get_time_hms([times[int(x)+ hor_cut]]), alt_range[int(y)+ ver_cut]))

# Click event which recognizes the scroll click to register the x, y values.
def onclick2(event):
    if event.button == 2:
       x, y = event.xdata, event.ydata
       try:
          #print(get_time_hms([times[int(x)]]), alt[int(y)])
          print(get_time_hms([x]))
       except Exception as e:
          print(e)
          return
       coordinates.append(np.where(times == x))
       time_alt.append((get_time_hms([x]), np.nan))



### Plots the reflectivity image
  # Parameters: 
  # refl        : Reflectivity matrix
  # t           : Time array 
  # fl_path     : Flight path array 
  # altrange    : Altitude at each cell
  # threshold   : Minimum value of the reflectivity
  # h_cut       : Cuts the image length (horizontal_cut)
  # v_cut       : Cuts the image height (used it to hide the reflectivities at the terrain-level to better see the colors in the area of interest)
  # x_step      : Step size between each x_tick on the x_axis. Reduce to see more times on the x_axis. 
  # y_step      : Step size between each y_tick on the y_axis. 
def plotImage(refl, t, fl_path, altrange, threshold, h_cut, v_cut, x_step, y_step): # 3000, 80, 120, 50
    ### Color map configuration
    #cmap = cm.get_cmap('gist_ncar_r') #gist_ncar_r, gist_ncar , and  prism are other colormaps 
    #cmap = gc("../colortables/myBkBlAqGrYeOrReViWh200.gp")
    #cmap = gc("../colortables/BlGrYeOrReVi200.gp")
    cmap = gc("../colortables/BlAqGrYeOrReVi200.rgb")


    ### masking the reflectivity image based on the given threshold. 
    masked_refl = refl # This is causing the rest of the plots to also have masked values according to the threshold.
    masked_refl[np.where(masked_refl < threshold)] = np.nan
    #masked_refl[np.where(masked_refl > 4)] = np.nan
    
    ### Plotting 
    fig = plt.figure(1)
    plt.imshow(masked_refl[h_cut:, v_cut:].transpose(), cmap = cmap, origin = 'lower')
    plt.colorbar()
    #plt.plot(np.array(fl_path[h_cut:])-v_cut) 
    plt.xticks(np.arange(0, len(t[h_cut:]), x_step), get_time_hms(t[np.arange(h_cut, len(t), x_step)]))
    plt.yticks(np.arange(0, len(altrange[v_cut:]), y_step), altrange[np.arange(v_cut, len(altrange), y_step)])
    plt.xlabel("Time (UTC)") 
    plt.ylabel("Altitude (m MSL)")
    plt.title(filename) 
    plt.grid(True)

    
    # Establishing connection between the figure and the onClick event. 
    cid = fig.canvas.mpl_connect('button_press_event', onclick)


### Plots the reflectivity vs time at the nearest gates with some reflectivity
  # Parameters: 
  # refl      : Reflectivity matrix 
  # r_c       : Number of gates above and below the flight level being searched for existing refl value
  # t         : Time array 
  # alt       : Flight alt at each time frame
  # altrange  : Altitude at each cell 
  # threshold : The min refl value, all values below it are masked
  # x_step    : Step size between two x_ticks, decrease value to increase the number of x_ticks on the x axis
def plotNearestGateReflectivity(refl, r_c, t, alt, altrange, threshold, x_step, h_cut):
    upper_gate = []
    lower_gate = [] 

    for i in range(len(alt)):
       abs_diff = np.abs(altrange - alt[i]) 
       flight_level = np.argmin(abs_diff) 
       refl_all_gates = refl.data[i,:]
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
    
    #print(len(lower_gate)) 
    #print(len(upper_gate)) 
    #plotXY(times, upper_gate)
    #plotXY(times, lower_gate)
    
    upper_gate = np.array(upper_gate) 
    lower_gate = np.array(lower_gate)

    global refl_upper_gate
    refl_upper_gate = np.array(upper_gate[:])
    refl_upper_gate[np.where(np.isnan(refl_upper_gate))] = mask_fill_value 
    global refl_lower_gate
    refl_lower_gate = np.array(lower_gate[:])
    refl_lower_gate[np.where(np.isnan(refl_lower_gate))] = mask_fill_value 


    if(threshold is not None):  
      upper_gate[np.where(np.isnan(upper_gate))] = mask_fill_value 
      lower_gate[np.where(np.isnan(lower_gate))] = mask_fill_value 
    
      lower_gate[lower_gate < threshold] = np.nan
      upper_gate[upper_gate < threshold] = np.nan
     
    #print(np.where(~(np.isnan(upper_gate))))
     
    fig2 = plt.figure(2)
    plt.plot(t[h_cut:], upper_gate[h_cut:], label="Closest Gate above flight track") 
    plt.plot(t[h_cut:], lower_gate[h_cut:], label="Closest Gate below flight track")
    plt.legend()
    plt.xlabel("Time (UTC)") 
    plt.ylabel("Reflectivity (dBZ)")
    plt.xticks(t[np.arange(h_cut, len(t), x_step)], get_time_hms(t[np.arange(h_cut, len(t), x_step)]))
    plt.grid(True) 
    plt.title(filename) 
    cid2 = fig2.canvas.mpl_connect('button_press_event', onclick2)


def plotDopplerVelocity(vel, t, r_c, alt, altrange, threshold, x_step, h_cut):
    upper_gate = []
    lower_gate = []
    #flight_gate = []

    for i in range(len(alt)):
       abs_diff = np.abs(altrange - alt[i]) 
       flight_level = np.argmin(abs_diff) 
       refl_all_gates = vel.data[i,:]
       refl_all_gates[:(flight_level-r_c)] = np.nan 
       refl_all_gates[(flight_level+r_c):] = np.nan
       refl_all_gates[refl_all_gates == mask_fill_value] = np.nan
       non_nan_refl_index = np.array(np.where(~(np.isnan(refl_all_gates)))).flatten()
       diff_gates = non_nan_refl_index - flight_level
       try: 
           upper_gate.append(vel[i][non_nan_refl_index[np.where(diff_gates == (np.min(diff_gates[diff_gates > 0])))[0][0]]])
       except Exception as e:
           #print("upper_gate: ", e) 
           upper_gate.append(np.nan) 
       try: 
           lower_gate.append(vel[i][non_nan_refl_index[np.where(diff_gates == (np.max(diff_gates[diff_gates < 0])))[0][0]]])
       except Exception as e:
           #print("lower_gate: ", e) 
           lower_gate.append(np.nan) 
       #try: 
       #    flight_gate.append(vel[i][flight_level])
       #except Exception as e:
       #    #print("flight_gate: ", e) 
       #    flight_gate.append(np.nan) 

    #print(len(lower_gate)) 
    #print(len(upper_gate)) 
    #plotXY(times, upper_gate)
    #plotXY(times, lower_gate)
    
    upper_gate = np.array(upper_gate) 
    lower_gate = np.array(lower_gate) 

    #flight_gate = np.array(flight_gate)
    if(threshold is not None):  
      lower_gate[refl_lower_gate < threshold] = np.nan
      upper_gate[refl_upper_gate < threshold] = np.nan
      #flight_gate[flight_gate < threshold] = np.nan
     
    fig3 = plt.figure(3)
    plt.plot(t[h_cut:], upper_gate[h_cut:], label="Closest Gate above flight track") 
    plt.plot(t[h_cut:], lower_gate[h_cut:], label="Closest Gate below flight track")
    #plt.plot(t, flight_gate, label = "UWKA fit measure")
    plt.legend()
    plt.xlabel("Time (UTC)") 
    plt.ylabel("Doppler Velocity (m/s)")
    plt.xticks(t[np.arange(h_cut, len(t), x_step)], get_time_hms(t[np.arange(h_cut, len(t), x_step)]))
    plt.grid(True) 
    plt.title(filename) 
    cid3 = fig3.canvas.mpl_connect('button_press_event', onclick2)


def plotLWC_Eddy(fil, t, h_cut, x_step):
    try:
        # Opens the netCdf file and exits if incorrect filename  is provided. 
        nc = netCDF4.Dataset(fil,'r')
    except:
        print("Please check the validity of the input files")
        quit()
    

    full_time = nc.variables['time'][:]
    dt = datetime(2017, 1, 1) 
    #print(str(dt)) 
    dt2 = datetime(1970, 1, 1) 
    #print(str(dt2))
    #print(str((dt - dt2).total_seconds())) 
    #print(str(time.mktime(dt.timetuple())))
    #print(dt.strftime("%s"))
    #print(t[len(t)-1])
    #print(full_time) 
    full_time += (dt - dt2).total_seconds() 
    full_time = np.array(full_time).astype(int)
    #print(full_time.shape) 
    st_index = np.where(full_time == int(t[0]))[0][0]
    #print(st_index)
    #quit()
    end_index = int(np.where(full_time == int(t[len(t)-1]))[0])+1
    #print(end_index)
    rosemount = nc.variables['rlwc'][st_index:end_index]
    rosemount[rosemount<0] = np.nan
    dmt100 = nc.variables['lwc100'][st_index:end_index]
    dmt100[dmt100 < 0] = np.nan
    dmtcdp = nc.variables['cdplwc_NRB'][st_index:end_index]
    dmtcdp[dmtcdp < 0] = np.nan
    nevz = nc.variables['nevlwc'][st_index:end_index]
    nevz[nevz < 0] = np.nan
    total_nevz = nc.variables['nevtwc'][st_index:end_index]
    total_nevz[total_nevz < 0] = np.nan
    eddy = nc.variables['turb'][st_index:end_index] 
    #eddy[eddy<0] = np.nan
    #print(np.array(rosemount).shape) 
    #print(np.array(t).shape)
    time_cut = np.where(full_time == int(t[h_cut]))[0][0] - st_index
    #print(time_cut) 
    fig4 = plt.figure(4) 
    time_plot = full_time[st_index:end_index]
    plt.plot(time_plot[time_cut:] , rosemount[time_cut:], label="Rosemount Icing Detector")
    plt.plot(time_plot[time_cut:], dmt100[time_cut:], label = "DMT100") 
    plt.plot(time_plot[time_cut:], dmtcdp[time_cut:], label = "DMT CDP") 
    plt.plot(time_plot[time_cut:], nevz[time_cut:], label = "Nevzorov") 
    plt.plot(time_plot[time_cut:], total_nevz[time_cut:], label = "Nev. Total LW") 
    plt.legend() 
    plt.xlabel("Time (UTC)") 
    plt.ylabel("LWC (g m^-3)") 
    plt.xticks(time_plot[np.arange(time_cut, len(time_plot), x_step)], get_time_hms(time_plot[np.arange(time_cut, len(time_plot), x_step)]))
    plt.grid(True) 
    plt.title(fil) 
    cid4 = fig4.canvas.mpl_connect('button_press_event', onclick2)
    
    fig5 = plt.figure(5) 
    plt.plot(time_plot[time_cut:], eddy[time_cut:]) 
    plt.ylabel("Eddy Dissipation Rate") 
    plt.xlabel("Time (UTC)") 
    plt.xticks(time_plot[np.arange(time_cut, len(time_plot), x_step)], get_time_hms(time_plot[np.arange(time_cut, len(time_plot), x_step)]))
    plt.grid(True)
    plt.title(fil) 
    cid5 = fig5.canvas.mpl_connect('button_press_event', onclick2)
   
    



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
    altitude = nc.variables['ALT'][:]
    alt_range = np.array(nc.variables['altrange'][:])
    vel_corr = nc.variables['velocity_corrected'][:]
    #print(vel_corr.shape) 

    ### Filling the flight_path array (global)
    for i in range(len(altitude)):
       abs_diff = np.abs(alt_range - altitude[i]) 
       flight_path.append(np.argmin(abs_diff))

    # plot refl image
    plotImage(reflectivity, times, flight_path, alt_range, min_refl, hor_cut, 80, 120, 50)

    # plot vel image 
    #plotImage(vel_corr, times, flight_path, alt_range, min_vel, hor_cut, 90, 120, 50) 

    # plot reflectivities of closest upper and lower gate
    plotNearestGateReflectivity(reflectivity, range_cut, times, altitude, alt_range, min_refl, 50, hor_cut)  

    plotDopplerVelocity(vel_corr, times, range_cut, altitude, alt_range, min_refl, 50, hor_cut)
    
    plotLWC_Eddy(argv[2], times, hor_cut, 50) 

   # Show all plots   
    plt.show() 

    # Close all connections
    fig.canvas.mpl_disconnect(cid)
    fig2.canvas.mpl_disconnect(cid2) 
    fig3.canvas.mpl_disconnect(cid3) 
    fig4.canvas.mpl_disconnect(cid4) 
    fig5.canvas.mpl_disconnect(cid5) 
