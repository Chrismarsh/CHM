#!/bin/bash
# This is needed becase crontab does not have same env variables are user
PATH=$PATH:/home/nwayand/custom/anaconda2:/home/nwayand/custom/anaconda2/bin:/home/nwayand/custom/anaconda2/bin

# Timer
start=$SECONDS

## Script paths
# Where scripts are located
ex_dir=/home/nwayand/python/HRDPS/
# CHM run dir
CHM_dir=/home/nwayand/snow_models/output_CHM/SnowCast/
# Config file to use
Configfile=Config_GEM_west_fortress.py

# Download GEM forecast
echo Downloading GEM
/home/nwayand/custom/anaconda2/bin/python $ex_dir"Download_HRDPS_GRIB2.py"  $ex_dir$Configfile

# Format grib2 to netcdf
echo Formating grib2 to netcdf
/home/nwayand/custom/anaconda2/bin/python $ex_dir"GRIB2_to_Netcdf.py" $ex_dir$Configfile

# Convert archived netcdf to CHM forcing
echo Converting archived netcdf to CHM ascii forcing files
/home/nwayand/custom/anaconda2/bin/python $ex_dir"Netcdf_to_CHM_forcing.py" $ex_dir$Configfile

# Run CHM for available forcing period
echo Running CHM
$CHM_dir"Run_CHM_all_GEM_hours.sh"

duration=$(( SECONDS - start ))
echo Took $duration seconds

