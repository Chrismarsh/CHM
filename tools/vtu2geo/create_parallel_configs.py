import vtk
from osgeo import gdal,ogr,osr
import math
import imp
import sys
import xml.etree.ElementTree as ET
import subprocess
import numpy as np
import os
import shutil
gdal.UseExceptions()  # Enable errors

def main():
    gdal.UseExceptions()  # Enable errors
    
    #####  load user configurable paramters here    #######
    # Check user defined configuraiton file
    if len(sys.argv) == 1:
        print('ERROR: main.py requires one argument [configuration file] (i.e. python main.py vtu2geo_config.py)')
        return

    # Get name of configuration file/module
    configfile = sys.argv[-1]

    # Load in configuration file as module
    X = imp.load_source('',configfile)
    
    # Grab variables
    input_path = X.input_path

    variables = X.variables
    parameters = []
    if hasattr(X,'parameters'):
        parameters = X.parameters

    #config_dir = 'parallel_configs'
    if hasattr(X,'config_dir'):
	config_dir = X.config_dir
	
    # Create folder if it doesn't exist
    if not os.path.exists(config_dir):
    	os.makedirs(config_dir)

    # Check if we want to constrain output to a example geotif
    constrain_flag = False
    if hasattr(X,'constrain_tif_file'):
        constrain_tif_file = X.constrain_tif_file
        var_resample_method = X.var_resample_method
        param_resample_method = X.param_resample_method
        constrain_flag = True

    output_path = input_path[:input_path.rfind('/')+1]
    if hasattr(X,'output_path'):
        output_path = X.output_path

    pixel_size = 0
    if hasattr(X,'pixel_size'):
        pixel_size = X.pixel_size

    user_define_extent = False
    if hasattr(X,'user_define_extent'):
        user_define_extent = X.user_define_extent

    # Get size for first rasterization (less than triangle min area)
    pixel_size = X.pixel_size

    #####
    reader = vtk.vtkXMLUnstructuredGridReader()
    pvd = ET.parse(input_path)
    pvd = pvd.findall(".//*[@file]")

    # Create execution file
    f2 = open(config_dir+'/parallel_commands.txt','w')

    # TODO: hardcoded
    para_ex = '/home/new365/CHM/tools/vtu2geo/main_parallel.py'

    for vtu in pvd:

	# Get full file name
        vtu_file  = vtu.get('file')
        
	# Copy templet config file
	new_config = config_dir+'/cfg_'+vtu_file+'.txt'
	shutil.copy(configfile, new_config)	

	# Open config file for writing
	f = open(new_config, 'a')

	# Write file name
	f.write('vtu_file="'+vtu_file+'"\n')

	# Close config file
	f.close()

	# add to master command file
        f2.write('python '+para_ex+' '+new_config+'\n')
    
    f2.close()

if __name__ == "__main__":

    main()
