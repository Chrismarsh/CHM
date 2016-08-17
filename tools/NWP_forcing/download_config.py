# Configuration file for Download_HRDPS_GRIB2.py

# Dir to put GEM grib2 files
download_dir = '/home/nwayand/data_sets/HRDPS/west'

# Forecast initializatio time (UTM)
# Options are: '00','06','12','18'
Init_H   = '00'

# Define HRDPS variables to download
Variable   = ['TMP_ISBL_1015','TMP_ISBL_1000','TMP_ISBL_0985','TMP_ISBL_0970','TMP_ISBL_0950','TMP_ISBL_0925',
'TMP_ISBL_0900','HGT_ISBL_1015','HGT_ISBL_1000','HGT_ISBL_0985','HGT_ISBL_0970','HGT_ISBL_0950','HGT_ISBL_0925','HGT_ISBL_0900',
'HGT_SFC_0','TMP_TGL_2','RH_TGL_2','WIND_TGL_10','WDIR_TGL_10','PRES_SFC_0','DLWRF_SFC_0','DSWRF_SFC_0','PRATE_SFC_0']



