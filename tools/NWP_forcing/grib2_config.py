#Configuration file for GRIB2_to_Netcdf.py

# Dir of grib2 files
download_dir = '/media/data2/GEM/west/grib2_current'

# Dir where output netcdf files go
output_dir   = '/media/data2/GEM/west/netcdf_archive'

# Dictionary list to change GEM variables names into CHM required names
var_dic = {'time':'datetime','TMP_P0_L103_GST0':'t','RH_P0_L103_GST0':'rh','WDIR_P0_L103_GST0':'vw_dir','WIND_P0_L103_GST0':'u','DLWRF_P8_L1_GST0_acc':'Qli','DSWRF_P8_L1_GST0_acc':'Qsi','PRATE_P0_L1_GST0':'p','PRES_P0_L1_GST0':'press','RPRATE_P8_L1_GST0_acc':'p_rain','SPRATE_P8_L1_GST0_acc':'p_snow','FPRATE_P8_L1_GST0_acc':'p_frz_rain','IPRATE_P8_L1_GST0_acc':'p_ice_pellet','ALBDO_P0_L1_GST0':'albedo','LHTFL_P0_L1_GST0':'latentHeat','SHTFL_P0_L1_GST0':'sensibleHeat'}

