import numpy as np 
import matplotlib
from matplotlib.colors import ListedColormap as lc 
from sys import argv

def getColorMap(filepath): 
    f = open(filepath, "r") 
    lines = f.readlines()
    found_rgb = False
    rgb = []
    for line in lines: 
        if found_rgb:
            values = line.split()
            values = map(int, values)
            #print(values) 
            rgb.append(values)
        #print(line.split())
        if (line.split() == ['#', 'r', 'g', 'b']):
            found_rgb = True
            #print(line)
    rgb = np.array(rgb)/255.0 
    return lc(rgb)



n_args = len(argv)
if n_args <= 1:
    print("please provide the filename(s) to be parsed") 
    quit()
else:
    filename = argv[1]
    try:
        # Opens colormap 
        rgb_values = getColorMap(filename) 
    except Exception as e:
        print(e) 
        #print("Please check the validity of the input files")
        quit()

#print(rgb_values) 

