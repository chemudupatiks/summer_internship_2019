import netCDF4
import numpy 
from sys import argv 


n_args = len(argv)
if n_args <= 1:
    print("please provide a filename") 
    quit()
else:
    filenames = argv[1:]
    try:
        for i in filenames:
            print("Variables in "+ i+ ": ")
            nc = netCDF4.Dataset(i, 'r')
            print(nc.variables.keys()) 
    except:
        print("Please check the validity of the input files")
        quit()

