# Configuration file for vtu2geo tool

# Vtu base name (i.e. base_name_XXX.vtu)
base = 'SC'

# pvd file path
input_path = '/Users/chris/Documents/PhD/code/CHM/output/fortress1.pvd'
output_path = '/home/chris/Documents/PhD/code/CHM/output_tif/'

# Output projection
EPSG=4326 # http://spatialreference.org/ref/epsg/
is_geographic = True
# Output variables
variables = ['t']  #set to None to dump all variables

# Output parameters
parameters = [] # paramters are one offs we want to extract from the vtu files

# Output pixel size that the mesh is interpolated to (?)
pixel_size = 300 # (m)

