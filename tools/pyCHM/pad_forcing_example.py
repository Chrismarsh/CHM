import os
import forcing_functions as ff

# Settings
agg_dt = 'H' # Output time step (stamp is beg of period)
input_dir    = os.path.normpath(r'C:\Users\new356\Model_Output\Forcing\Pad_test\input')
output_dir   = os.path.normpath(r'C:\Users\new356\Model_Output\Forcing\Pad_test\output')

# Pad forcing files
ff.pad_forcing_to_common_length(agg_dt,input_dir,output_dir)

