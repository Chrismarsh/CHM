# Configuration file for vtu2geo tool

# pvd file path
input_path = '/Users/chris/Documents/PhD/code/CHM/output/meshes/marmot_deform.pvd'

#defaults to input_path if not specified
#output_path = '/Users/chris/Documents/PhD/code/CHM/output_tif'

# Output variables
variables = ['t']  #set to None to dump all variables

# Output parameters
parameters = ['Elevation'] # parameters are one offs we want to extract from the vtu files, only will be written out once

# Output pixel size that the mesh is interpolated to
# defaults to (min+max)/2 if not specified
pixel_size = 10 # (m)

#### Options for outputing second clipped geotiffs ####

# (Optional) file to constrain 2nd output tif. Comment out if not used.
# Uses resolutoin from this tif
# Uses extent from this tif, unless specified below with user_define_extent flag
constrain_tif_file = '/home/nwayand/data_sets/Remote/Landsat/LandsatNDSI/NDSI140928A_pro.tif'

## Define output info
user_define_extent = True # If True, below extent is used to clip,
# otherwise extent is taken from the constrain_tif_file extent
o_xmin = -1287634.12811
o_xmax = -1269158.43038
o_ymin = 1394814.5179
o_ymax = 1417759.53926

# Resample methods for second geotiff (one for each variable/parameter)
var_resample_method    = {'t':'average'} # Methods to use when calculating clipped raster
param_resample_method  = {'Elevation':'average'} # Methods to use when calculating clipped raster


