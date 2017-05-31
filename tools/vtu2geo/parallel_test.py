# Configuration file for vtu2geo tool

# Vtu base name (i.e. base_name_XXX.vtu)
base = 'SC'

# pvd file path
input_path = '/home/new365/config_CHM/BlowingAvaSnow/OBS_no_blow_yes_ava_filled_dist/meshes/SC.pvd'
output_path = '/home/new365/config_CHM/BlowingAvaSnow/OBS_no_blow_yes_ava_filled_dist/meshes/'
config_dir = '/home/new365/CHM/tools/vtu2geo/para_configs'

constrain_tif_file = '/home/new365/config_CHM/modscag_crho.tif'

# Pixel size for first raserization (must be smaller than min_area of triangles)
pixel_size = 5 # m

## Define output info
user_define_extent = True # If True, below extent is used to clip,
# otherwise extetn is taken from the constrain_tif_file extent
# Extent (Fortress domain)
o_xmin = -1287634.12811
o_xmax = -1269158.43038
o_ymin = 1394814.5179
o_ymax = 1417759.53926
# # Extent (CRHO domain)
# o_xmin = -1374651.2329267
# o_xmax = -1225383.09324314
# o_ymin = 1381800.81154171
# o_ymax = 1537266.50438795

# Output variables
variables = ['snowcoverfraction']  #set to None to dump all variables
var_resample_method    = {'snowcoverfraction':'average'} # Methods to use when calculating clipped raster

# Output parameters
parameters = ['Elevation'] # paramters are one offs we want to extract from the vtu files
param_resample_method    = {'[param] landcover':'mode','Elevation':'average'} # Methods to use when calculating clipped raster

