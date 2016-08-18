# Configuration file for vtu2geo tool

# Vtu base name (i.e. base_name_XXX.vtu)
base = 'SC'

# Input path to where output vtu files are located
input_path = '/home/nwayand/snow_models/output_CHM/SnowCast/output/'
# output_path = '/home/chris/Documents/PhD/code/CHM/output_tif/'

# Output projection
EPSG=26911 # http://spatialreference.org/ref/epsg/

# Output variables
variables = ['t','p']  #set to None to dump all variables

# Output parameters
parameters = [] # paramters are one offs we want to extract from the vtu files

# Output pixel size that the mesh is interpolated to (?)
pixel_size = 10 # (m)

