import xarray as xr
import os
import glob
import imp
import sys
import numpy as np
import pandas as pd
import datetime
import json
import time
import utm
start_time = time.time()
# Hack to force datetimes to display in GMT/UTC (numpy 1.11.1 has fixed this but other dependent modules (pynio) can't handel numpy 1.11.1)
os.environ['TZ'] = 'GMT'
time.tzset()

# Load in config file
#######  load user configurable paramters here    #######
# Check user defined configuraiton file
if len(sys.argv) == 1:
    sys.error('GRIB2_to_CHM_forcing.py requires one argument [configuration file] (i.e. python GRIB2_to_CHM_forcing.py forcing_config.py')


# Get name of configuration file/module
configfile = sys.argv[-1]

# Load in configuration file as module
X = imp.load_source('',configfile)

# Assinge to local variables
netcdf_dir = X.netcdf_dir
ascii_dir   = X.ascii_dir
Forcing_config_file = X.Forcing_config_file
lat_r = X.lat_r
lon_r = X.lon_r
local_time_offset = X.local_time_offset
coordsystem = X.coordsystem
if coordsystem=='pro':
    utm_zone = X.utm_zone

# Move to input
os.chdir(netcdf_dir)

# Get all file names
all_files =  glob.glob('*.nc')

# Sort file names (open_mfdataset _should_ handle this, but this is a hack to put in correct order)
all_files = sorted(all_files)

# Load in all netcdf files
def preprocess(x):
    # Get the end of the most recent forecast (in UTC)
    c_for_end = np.datetime64(datetime.datetime.strptime(time.strftime("%Y%m%d"),"%Y%m%d")+datetime.timedelta(days=2)) # So ugly, need to make cleaner call
    # Check if it is the most recent forecast
    # If so grab 48hr forecast
    if (x.datetime[-1]==c_for_end):
	x = x.isel(datetime=np.arange(1,47))
    # Otherwise only grab the 12 hours
    else:
    	x = x.isel(datetime=np.arange(1,25)) # Grab forecast hours 02 - 25(01) (we can't use the first 01 forecsat hour, because radiation vars were saved accumulated, thus we don't have the first value. This is fine we just use later forecast hours by 1 hour
    
    x.load()
    return x

ds = xr.open_mfdataset(all_files,concat_dim='datetime',engine='netcdf4',preprocess=lambda x: preprocess(x))

# Adjust to local time zone (i.e. from UTC to MST, local_time_offset should = -7)
ds['datetime'] = pd.to_datetime(ds.datetime.values) + datetime.timedelta(hours=local_time_offset)

# Move to ascii dir
if not os.path.isdir(ascii_dir):
    os.mkdir(ascii_dir)
os.chdir(ascii_dir)

# Extract grid cells we want to export
print 'Extracting cells within lat/long box'
ds = ds.where((ds.gridlat_0>lat_r[0]) & (ds.gridlat_0<lat_r[1]) & (ds.gridlon_0>lon_r[0]) & (ds.gridlon_0<lon_r[1]), drop=True)

# Get small dataset of lat and long
grid = ds[['gridlat_0', 'gridlon_0','HGT_P0_L1_GST0']]

# Initalize dictionary of metadata
meta = {}

# Loop through each point within our lat long box
for i in range(ds.coords['xgrid_0'].size):
    for j in range(ds.coords['ygrid_0'].size):
	sub_grid = grid.isel(xgrid_0=i, ygrid_0=j)
        # Is this test needed here???
	if ((sub_grid.gridlat_0>lat_r[0]) & (sub_grid.gridlat_0<lat_r[1]) & (sub_grid.gridlon_0>lon_r[0]) & (sub_grid.gridlon_0<lon_r[1])):	
            # Select point
	    sub_ds = ds.isel(xgrid_0=i, ygrid_0=j)
            # Convert to dataframe
	    df = sub_ds.to_dataframe()
	    # Drop unwanted vars/coords 
            df.drop(['HGT_P0_L1_GST0','xgrid_0','ygrid_0','gridlat_0','gridlon_0'], axis=1, inplace=True)
	    # Write to ascii CHM format
 	    df.to_csv('point_'+str(i)+'_'+str(j)+'.chm',sep='\t',date_format='%Y%m%dT%H%M%S')
	    # Build dicting.json',cnary of point metadata
	    coords = []
	    if coordsystem=='geo':
		coords.append(float(sub_grid.gridlon_0.values))
	 	coords.append(float(sub_grid.gridlat_0.values))
	    elif coordsystem=='pro':
	    	coords = utm.from_latlon(sub_grid.gridlat_0,sub_grid.gridlon_0, utm_zone)
            else:
		print 'Unknown coords = ' + coordsystem

	    meta['point_'+str(i)+'_'+str(j)] = {'file':ascii_dir+'/point_'+str(i)+'_'+str(j)+'.chm',
		"longitude":coords[0],
		"latitude":coords[1],
		"elevation":sub_grid.HGT_P0_L1_GST0.values[0],
		"filter": {
        	"scale_wind_speed": {
          	"variable": "u",
          	"Z_F": 10
        	}}}

# Write out json file with metadata for these points
print("Writing out CHM forcing json file")
with open('GEM_forcing.json', 'w') as outfile:
    json.dump(meta, outfile,indent=4, sort_keys=True)

print("--- %s minutes ---" % ((time.time() - start_time)/60))


print 'finished!'
