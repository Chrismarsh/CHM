################################################
# Paths
################################################
# Dir to put GEM grib2 files
download_dir = '/media/data2/GEM/west/grib2_current'

# Dir where output netcdf files go
netcdf_dir  = '/media/data2/GEM/west/netcdf_archive'

# Dir where output ascii files should go
ascii_dir    = '/media/data2/GEM/west/CHM_archive_fortress'

################################################
# Configuration for Download_HRDPS_GRIB2.py
################################################

# Forecast initializatio time (UTM)
# Options are: '00','06','12','18'
Init_H   = '00'

domain = 'west'

# Define HRDPS variables to download (names match file names in HRDPS system)
Variable   = ['TMP_ISBL_1015','TMP_ISBL_1000','TMP_ISBL_0985','TMP_ISBL_0970','TMP_ISBL_0950','TMP_ISBL_0925',
'TMP_ISBL_0900','HGT_ISBL_1015','HGT_ISBL_1000','HGT_ISBL_0985','HGT_ISBL_0970','HGT_ISBL_0950','HGT_ISBL_0925','HGT_ISBL_0900',
'HGT_SFC_0','TMP_TGL_2','RH_TGL_2','WIND_TGL_10','WDIR_TGL_10','PRES_SFC_0','DLWRF_SFC_0','DSWRF_SFC_0','PRATE_SFC_0',
'WEASN_SFC_0','WEARN_SFC_0','WEAPE_SFC_0','WEAFR_SFC_0','LHTFL_SFC_0','SHTFL_SFC_0','ALBDO_SFC_0']

##########################################
#Configuration for GRIB2_to_Netcdf.py
##########################################

# Dictionary list to change GEM variables names into CHM required names (Names match internal variable names within grib files (note: they don't match above Variable names!))
var_dic = {'time':'datetime','TMP_P0_L103_GST0':'t','RH_P0_L103_GST0':'rh','WDIR_P0_L103_GST0':'vw_dir','WIND_P0_L103_GST0':'u','DLWRF_P8_L1_GST0_acc':'Qli','DSWRF_P8_L1_GST0_acc':'Qsi','PRATE_P0_L1_GST0':'p','PRES_P0_L1_GST0':'press','RPRATE_P8_L1_GST0_acc':'p_rain','SPRATE_P8_L1_GST0_acc':'p_snow','FPRATE_P8_L1_GST0_acc':'p_frz_rain','IPRATE_P8_L1_GST0_acc':'p_ice_pellet','ALBDO_P0_L1_GST0':'albedo','LHTFL_P0_L1_GST0':'latentHeat','SHTFL_P0_L1_GST0':'sensibleHeat'}

#########################################
# Configuration for Netcdf_to_CHM.py
#########################################

# Name for output .json file containing metadata for forcing files
Forcing_config_file = 'GEM_forcing.json'

# Offset from UTM to local time (i.e. Mountain standard time = -7)
# CHM forcing files will be in this time zone
local_time_offset = 0

# Coordinate system for output CHM files
coordsystem = 'geo' # 'geo' = geographic, 'pro' = projected
utm_zone = 11 # Only used if coordsystem='pro'

#### [min,max] extents for bounding box of latitude and longitude

# Fortress
lat_r = [50.79,50.876]
lon_r = [-115.26,115.151]

