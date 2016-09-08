# Configuration file for Netcdf_to_CHM.py

# Dir of netcdf files
netcdf_dir = '/media/data2/GEM/west/netcdf_archive'

# Dir where output ascii/netcdf/json files should go
output_dir   = '/media/data2/GEM/west/CHM_archive'

# Name for output .json file containing metadata for forcing files
Forcing_config_file = 'GEM_forcing.json'

# Offset from UTM to local time (i.e. Mountain standard time = -7)
# CHM forcing files will be in this time zone
local_time_offset = -7

# Coordinate system for output CHM files
coordsystem = 'geo' # 'geo' = geographic, 'pro' = projected
utm_zone = 11 # Only used if coordsystem='pro'

#### [min,max] extents for bounding box of latitude and longitude

# Bow river basin
#lat_r = [50.411581,51.218712]
#lon_r = [-115.793152,-114.362183]

# Entire GEM west domain
lat_r = [0,90]
lon_r = [-180,180]


