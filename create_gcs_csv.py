import pandas as pd
import numpy as np 
import netCDF4 
import time 
from datetime import datetime

def get_consecutive_num(array):
    cons = []
    sub = []
    for i in array:
        if len(sub) == 0:
            sub.append(i)
        elif i-sub[-1] == 1:
            sub.append(i) 
        else:
            cons.append(sub) 
            sub = [i]
    cons.append(sub) 
    return cons
        



# above : True or False. True implies a csv is generated based on data above the filght level, else data below 
def generate_gcs_csv(above):
    gcs_csv_all = [] 
    columns = ['filename', 'iop', 'date', 'group_num', 'cell_num','start_index', 'end_index', 'start_seconds', 'end_seconds', 'start_time', 'end_time', 'threshold']
    group_csv = pd.read_csv("../csv/groups.csv")
    #print(group_csv.columns)
    mask_fill_value = -32767
    main_groups = group_csv.loc[group_csv['main'] == 1]
    #print(group_csv) 
    #print(main_groups)
    for j in range(main_groups.shape[0]):
        #print(main_groups['filename'].iloc[j]) 
        try: 
            nc = netCDF4.Dataset("../SNOWIE/" + main_groups['filename'].iloc[j]) 
        except Exception as e:
            print(e) 
            print("check the filenames in groups.csv") 
            quit()
        time_array = nc.variables['time'][:]
        altrange = nc.variables['altrange'][:]
        alt = nc.variables['ALT'][:]
        refl = nc.variables['reflectivity'][:]
        date = main_groups['date'].iloc[j]
        filename = main_groups['filename'].iloc[j]
        iop = main_groups['iop'].iloc[j]
        group_num = main_groups['group_num'].iloc[j]
        dt = datetime(1970, 1, 1) 
        start_time = date +" "+ main_groups['start_time'].iloc[j]
        end_time = date +" "+ main_groups['end_time'].iloc[j]
        start_time = (datetime.fromtimestamp(time.mktime(time.strptime(start_time, "%m/%d/%Y %H.%M.%S")))-dt).total_seconds()
        end_time = (datetime.fromtimestamp(time.mktime(time.strptime(end_time, "%m/%d/%Y %H.%M.%S")))-dt).total_seconds()
        st_index = np.where(time_array.astype(int) == start_time)[0][0]
        end_index = np.where(time_array.astype(int) == end_time)[0][-1]
        #print(st_index)
        #print(start_time)
        #print(end_index) 
        #print(end_time)
        #print(date) 
        r_c = 5
        threshold = main_groups['threshold'].iloc[j]
        #print(get_time_intervals(time, main_groups['threshold'].iloc[i])) 
        #print(nc.variables)
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
        upper_gate[np.where(np.isnan(upper_gate))] = mask_fill_value 
        lower_gate[np.where(np.isnan(lower_gate))] = mask_fill_value
        #upper_gate[upper_gate < threshold] = np.nan
        #lower_gate[lower_gate < threshold] = np.nan
        upper_gate = np.where(upper_gate >= threshold)
        lower_gate = np.where(lower_gate >= threshold)
        if above:
            index_in_range = np.take(upper_gate, np.where(np.logical_and(upper_gate >= st_index, upper_gate <= end_index))[1])
        else: index_in_range  = np.take(lower_gate, np.where(np.logical_and(lower_gate >= st_index, lower_gate <= end_index))[1])
        #print(index_in_range) 
        gcs_csv = []
        cons = get_consecutive_num(index_in_range) 
        for k in range(len(cons)):
            st_time = time_array[cons[k][0]]
            e_time = time_array[cons[k][-1]]
            gcs_csv.append([filename, iop, date, group_num, k+1, cons[k][0], cons[k][-1], st_time, e_time, time.strftime('%H:%M:%S', time.gmtime(st_time)), time.strftime('%H:%M:%S', time.gmtime(e_time)), threshold])
            #print(time.strftime('%H:%M:%S', time.gmtime(st_time)))
            #print(st_time) 
            #print(time.strftime('%H:%M:%S', time.gmtime(e_time)))
            #print(e_time)
        #print(gcs_csv) 
        #print(upper_gate) 
        #print(lower_gate)
        gcs_csv_all+= gcs_csv 
    gcs_csv_all = pd.DataFrame(gcs_csv_all, columns = columns) 
    return gcs_csv_all


gcs = generate_gcs_csv(False)
#print(gcs)
try:
    msg = gcs.to_csv("../csv/gcs.csv") 
except Exception as e:
    print(e) 
    quit()
print("gcs.csv created in csv directory!") 
#print(msg)
