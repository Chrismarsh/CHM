#!/home/nwayand/custom/anaconda2/bin/python

import xarray as xr
import os
import imp
import sys
import numpy as np
import pandas as pd
import datetime
import json
import time
import utm
import threading
start_time = time.time()

def add_datetime_dataset(ds_in,local_time_offset):
	import re
        # Replace dummey time variable with datetime64
        init_time = ds_in.RH_P0_L103_GST0.attrs['initial_time']
        DT = re.split(' ',init_time)
        D = re.split('/',DT[0])
        HHMM = re.split(':',DT[1].replace("(","").replace(")",""))
        T_init = datetime.datetime(int(D[2]),int(D[0]),int(D[1]),int(HHMM[0])+1-local_time_offset,int(HHMM[1])) # +1 is because we start on forecast houor 1 (disregard hour 00 as it is just initial conditions... and no precip values)
        TimeS = pd.date_range(T_init, periods=len(ds_in.forecast_hour),freq='H')
        ds_in['forecast_hour'] = TimeS
	ds_in.rename({'forecast_hour':'time'},inplace=True)
        return ds_in

def add_datetime_dataarray(dr_in,local_time_offset):
        import re
        # Replace dummey time variable with datetime64
        init_time = dr_in.initial_time
        DT = re.split(' ',init_time)
        D = re.split('/',DT[0])
        HHMM = re.split(':',DT[1].replace("(","").replace(")",""))
        T_init = datetime.datetime(int(D[2]),int(D[0]),int(D[1]),int(HHMM[0])+1-local_time_offset,int(HHMM[1])) # +1 is because we start on forecast houor 1 (disregard hour 00 as it is just initial conditions... and no precip values)
        TimeS = pd.date_range(T_init, periods=len(dr_in.forecast_hour),freq='H')
        dr_in['forecast_hour'] = TimeS
	dr_in = dr_in.rename({'forecast_hour':'time'})
        return dr_in

# Create preprocessing function to handle loading and time
def preprocess(x):
    x.load()
    x = x.drop('gridrot_0')
    [cvar for cvar in x.data_vars]
    x = xr.concat([x],dim='forecast_hour')
    # units can be in hour or minuts (I have no idea why)
    if x[cvar].forecast_time_units=='hours':
    	x['forecast_hour'] = x[cvar].attrs['forecast_time']
    elif x[cvar].forecast_time_units=='minutes':
   	x['forecast_hour'] = x[cvar].attrs['forecast_time']/60 # minutes to hours
    else:
	import sys
	print('found forecast_time_units not expected')
	sys.exit()
    return x

def load_GEM_4d_var(PresLevs,UA_files,var_name,var_name_new,preprocess):
	import xarray as xr
        ds_UA = xr.Dataset()
        for cP in PresLevs:
                # Get all time step files for this pressure level
                cfiles = [x for x in UA_files if cP in x]
                #print(cfiles)
                # Load and conact by time
                #ds_t = xr.open_mfdataset(cfiles,concat_dim='time',engine='pynio',lock=threading.Lock())
                ds_t = xr.open_mfdataset(cfiles,concat_dim='forecast_hour',engine='pynio',preprocess=lambda x: preprocess(x))
		# Rename var based on file name pressure level
                ds_t.rename({var_name:var_name_new+cP},inplace=True)
                # Combine
                ds_UA = xr.merge([ds_UA,ds_t])
        # Drop unneed vars
        # ds_UA = ds_UA.drop('gridrot_0')
        # Conact by pressure level
        ds_com = xr.concat([ds_UA[cvar] for cvar in ds_UA.data_vars],dim='PressLev')
        # Rename
	ds_com.name = var_name_new
	# return
        return ds_com

# Start Main

# Load in config file
#######  load user configurable paramters here    #######
# Check user defined configuraiton file
if len(sys.argv) == 1:
    sys.error('GRIB2_to_CHM_forcing.py requires one argument [configuration file] (i.e. GRIB2_to_CHM_forcing.py forcing_config.py')

# Get name of configuration file/module
configfile = sys.argv[-1]

# Load in configuration file as module
X = imp.load_source('',configfile)

# Assinge to local variables
download_dir = X.download_dir
output_dir   = X.output_dir
local_time_offset = X.local_time_offset
lat_r = X.lat_r
lon_r = X.lon_r
export_netcdf = X.export_netcdf

# Move to input
os.chdir(download_dir)

# Get all file names
all_files = os.listdir(download_dir)

# import and combine all grib2 files
print 'Opening all grib2 data files'

# Option, create different file list for surface and upper air files. Load in separatly. 
UA_TMP = [x for x in all_files if 'TMP_ISBL' in x]
UA_HGT = [x for x in all_files if 'HGT_ISBL' in x]

# Presure levels to extract air temperature from
PresLevs = ['1015','1000','0985','0970','0950','0925','0900']

srf_files = [x for x in all_files if not 'HGT_ISBL' in x and not 'TMP_ISBL' in x]

# Load surface variables
print('Loading Surface variables')
#ds = xr.open_mfdataset(srf_files,concat_dim='time',engine='pynio',lock=threading.Lock())
ds = xr.open_mfdataset(srf_files,concat_dim='forecast_hour',engine='pynio',preprocess=lambda x: preprocess(x))
ds = add_datetime_dataset(ds,local_time_offset)

# Load upper atmosphere variables
print('Loading upper air Temperature')
ds_UA_T   = load_GEM_4d_var(PresLevs,UA_TMP,'TMP_P0_L100_GST0','TMP_',preprocess)
ds_UA_T   = add_datetime_dataarray(ds_UA_T,local_time_offset)

print('Loading upper air height at pressurelevels')
ds_UA_HGT = load_GEM_4d_var(PresLevs,UA_HGT,'HGT_P0_L100_GST0','HGT_',preprocess)
# Convert Geopotential height to geometric height (http://www.pdas.com/geopot.pdf)
#ds_UA_HGT = ds_UA_HGT* 6371*1000 / (6371*1000/()-ds_UA_HGT)
ds_UA_HGT = add_datetime_dataarray(ds_UA_HGT,local_time_offset)

# Merge together
ds_UA = xr.merge([ds_UA_T,ds_UA_HGT])

# Approx method of calculating lapse rate (diff from lower and upper atmos temp)
ds_LR = -1*(ds_UA.TMP_[6,:,:,:] - ds_UA.TMP_[0,:,:,:]) / (ds_UA.HGT_[6,:,:,:] - ds_UA.HGT_[0,:,:,:])
ds_LR.name = 't_lapse_rate'

# More accurate but far too slow method of calc lapse rate
#print("Calculating lower-atmosphere lapse rate")
#ds_LR = ds.TMP_P0_L103_GST0 * 0 - 9999 # Quick way to make dataarray with -9999 values but with correct dims/coords
#for cts in np.arange(0,len(ds_UA.time)):
#        print(cts)
#	for cx in ds_UA.xgrid_0.values:
#		for cy in ds_UA.ygrid_0.values:
#			x_temp = ds_UA.TMP_[:,cts,cy,cx] # Grab the vertical profile of air temperature
#			y_hgt  = ds_UA.HGT_[:,cts,cy,cx] # Grab the vertical heights of air temperature values
#			s      = np.polyfit(y_hgt,x_temp,1) # Fit a line to the data
#			ds_LR[cts,cy,cx].values = s[0] # Grab the slope (first element)

ds = xr.merge([ds,ds_LR])

# Get time step (in seconds)
dt_s = 1*60*60 # seconds (hard coded for now)

print 'Renaming variables to CHM syntax'
# Change names
var_dic = {'time':'datetime','TMP_P0_L103_GST0':'t','RH_P0_L103_GST0':'rh','WDIR_P0_L103_GST0':'vw_dir','WIND_P0_L103_GST0':'u','DLWRF_P8_L1_GST0_acc':'Qli','DSWRF_P8_L1_GST0_acc':'Qsi','PRATE_P0_L1_GST0':'p','PRES_P0_L1_GST0':'press'}
ds.rename(var_dic,inplace=True)

print 'Converting units to CHM units'
##### Convert units to CHM requirements

# Air temperature
ds['t'] = ds.t - 273.15

# Precipitation rate to accumulation (mm)
ds['p'] = ds['p'] * dt_s # density of water cancels out m to mm conversion

# Pressure
ds['press'] = ds['press'] / 100 # Pa to hPa

##### Radiation
# Shortwave radiation incoming
Qsi_wm2 = ds.Qsi.diff(dim='datetime')/dt_s # j/m2 to j/(s*m2)
# Set SW values just below zero to zero
Qsi_wm2.values[Qsi_wm2.values<0] = 0
# First value is unknown (downside of saving as accum...) so we set it to -9999
ds['Qsi'] = xr.concat([ds.Qsi[0,:,:]*0-9999,Qsi_wm2],dim='datetime').transpose('datetime','ygrid_0','xgrid_0')

# Longwave radiation incoming
Qli_wm2 = ds.Qli.diff(dim='datetime')/dt_s # j/m2 to j/(s*m2)
# First value is unknown (downside of saving as accum...) so we set it to -9999
ds['Qli'] = xr.concat([ds.Qli[0,:,:]*0-9999,Qli_wm2],dim='datetime').transpose('datetime','ygrid_0','xgrid_0')

# Move to output dir
os.chdir(output_dir)

# REMOVE ALL previous files
print('Deleting previous forcing files and config file!')
files = os.listdir(output_dir)
for file in files:
    os.remove(file)

# Export to netcdf (optional)
if export_netcdf:
    print('Writing netcdf file')
    ds.to_netcdf('GEM_2_5km_west.nc',engine='netcdf4')

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
	    utmcoords = utm.from_latlon(sub_grid.gridlat_0,sub_grid.gridlon_0, 11)
	    meta['point_'+str(i)+'_'+str(j)] = {'file':output_dir+'/point_'+str(i)+'_'+str(j)+'.chm',
		"easting":utmcoords[0],
		"northing":utmcoords[1],
		"elevation":sub_grid.HGT_P0_L1_GST0.values[0]}

# Write out json file with metadata for these points
print("Writing out CHM forcing json file")
with open('GEM_forcing.json', 'w') as outfile:
    json.dump(meta, outfile,indent=4, sort_keys=True)

print("--- %s minutes ---" % ((time.time() - start_time)/60))


print 'finished!'

