# Configuration file for vtu2geo tool

# pvd file path
input_path = '/Users/chris/Documents/PhD/code/CHM/output/meshes/marmot_deform.pvd'

#defaults to input_path if not specified
#output_path = '/Users/chris/Documents/PhD/code/CHM/output_tif'

# Output variables
variables = ['t']  #set to None to dump all variables

# Output parameters
parameters = [] # parameters are one offs we want to extract from the vtu files, only will be written out once

# Output pixel size that the mesh is interpolated to
# defaults to (min+max)/2 if not specified
pixel_size = 10 # (m)

