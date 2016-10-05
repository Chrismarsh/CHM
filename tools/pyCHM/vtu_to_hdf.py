import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import glob
import time
import vtk
from vtk.util import numpy_support as VN
import matplotlib.tri as tri

import vtu_functions as vfunc

# Dir to vtu files
vtu_dir   = os.path.normpath(r'/home/nwayand/snow_models/output_CHM/k_country/Default_run/meshes')  
prefix    = 'GemK'
hdf_dir   = os.path.normpath(r'/home/nwayand/snow_models/output_CHM/k_country/Default_run')

# Move to vtu dir
os.chdir(vtu_dir)
 
# Get list of files
vtu_files = glob.glob(prefix+'*.vtu') 

# Get time stamps
time_stamp = vfunc.get_vtu_time(vtu_files,prefix)

# Save triangle info
tri_file = hdf_dir + '/' + 'triangle_info.npy'
vfunc.save_triangle_info(vtu_files[0],tri_file)

# Save multple variables
all_vars = ['t','p']
vfunc.vtu_to_hdf(all_vars,vtu_files,time_stamp,hdf_dir)


