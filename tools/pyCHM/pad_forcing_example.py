import os
import forcing_functions as ff

# Settings
agg_dt = 'H' # Output time step (stamp is beg of period)
input_dir    = '/media/data2/nicway/forcing/input'
output_dir   = '/media/data2/nicway/forcing/output'

# Pad forcing files
ff.pad_forcing_to_common_length(agg_dt,input_dir,output_dir)

