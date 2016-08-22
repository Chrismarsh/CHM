#Configuration file for GRIB2_to_CHM_forcing.py

# Dir of grib2 files
download_dir = '/home/nwayand/data_sets/HRDPS/west'

# Dir where output ascii/netcdf/json files should go
output_dir   = '/home/nwayand/data_sets/HRDPS/forcing'

# Option to export to netcdf file (in addition to CHM formated ascii files)
export_netcdf = False

# Name for output .json file containing metadata for forcing files
Forcing_config_file = 'GEM_forcing.json'

# Offset from UTM to local time (i.e. Mountain standard time = -6)
local_time_offset = -6

# Dictionary list to change GEM variables names into CHM required names
var_dic = {'time':'datetime','TMP_P0_L103_GST0':'t','RH_P0_L103_GST0':'rh','WDIR_P0_L103_GST0':'vw_dir','WIND_P0_L103_GST0':'u','DLWRF_P8_L1_GST0_acc':'Qli','DSWRF_P8_L1_GST0_acc':'Qsi','PRATE_P0_L1_GST0':'p','PRES_P0_L1_GST0':'press'}

#### [min,max] extents for bounding box of latitude and longitude

# Bow river basin
lat_r = [50.411581,51.218712]
lon_r = [-115.793152,-114.362183]
# Fortress and Marmot
#lat_r = [50.720319,51.099587]
#lon_r = [-115.41481,-114.941025]
# Fortress
#lat_r = [50.737705,50.90687]
#lon_r = [-115.366745,-115.002823]

