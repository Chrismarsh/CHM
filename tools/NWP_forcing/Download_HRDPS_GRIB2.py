#!/home/nwayand/custom/anaconda3/bin/python

import os
import numpy as np
import wget
import sys
import imp
import time

# Load in config file
#######  load user configurable paramters here    #######
# Check user defined configuraiton file 
if len(sys.argv) == 1:
    sys.error('Download_HRDPS_GRIB2.py requires one argument [configuration file] (i.e. Download_HRDPS_GRIB2.py download_config.py')

# Get name of configuration file/module
configfile = sys.argv[-1]

# Load in configuration file as module
X = imp.load_source('',configfile)

# Assinge to local variables
download_dir = X.download_dir
Init_H = X.Init_H
Variable = X.Variable

# Move to input
os.chdir(download_dir)

# REMOVE ALL previous files
print('Removing ALL previous grib2 files')
files = os.listdir()
for file in files:
    os.remove(file)

# Define format of date in file names
YYYYMMDDHH = time.strftime("%Y%m%d")+Init_H

# Fixed parameters for path and file name syntax
base_url = 'http://dd.weather.gc.ca/model_hrdps'
filetype = 'grib2'
Forc_H   = np.arange(1,49)
sep      = '/'
office     = 'CMC'
system     = 'hrdps'
domain     = 'west'
projection = 'ps2.5km'
Phhh       = 'P'
mm         = '-00'
ending     = '.grib2'
us         = '_'

# For each forecast hour (Forc_H)
for c_Forc_H in Forc_H:
    # Convert to string and pad with zeros
    FH_s = np.array_str(c_Forc_H).zfill(3)

    # Get current path
    cpath = base_url + sep + domain + sep + filetype + sep + Init_H + sep + FH_s + sep
    print(cpath)

    # For each variable
    for cVar in np.arange(0,len(Variable)):
        
        # Make file path
        cfile = office + us + system + us + domain + us + Variable[cVar] + us + projection + us + YYYYMMDDHH + us + Phhh + FH_s +  mm + ending

        # download file
        filename = wget.download(cpath+cfile)

