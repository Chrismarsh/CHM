## Functions to work with CHM forcing files

def pad_forcing_to_common_length(agg_dt,input_dir,output_dir):
    import os
    import sys
    import glob
    import datetime
    import pandas as pd
    import numpy as np
    
    # Get input files
    os.chdir(input_dir)
    in_files = glob.glob('*')

    # Loop over files to get earliest/latest dates
    date_min = datetime.datetime(3000,1,1)
    date_max = datetime.datetime(1000,1,1)
    for cf in in_files:
        # load it in
        c_df = pd.read_csv(cf,sep='\t',parse_dates=True) 
        c_df.set_index('datetime',inplace=True)
        c_df.index = pd.to_datetime(c_df.index)
        # start/end dates
        c_min = c_df.index[0]
        c_max = c_df.index[-1]
        # keep track of dates
        if c_min < date_min:
            date_min = c_min
        if c_max > date_max:
            date_max = c_max

    print('Earliest date is: ' + str(date_min))
    print('Earliest date is: ' + str(date_max))

    new_index = pd.DatetimeIndex(start=date_min,end=date_max,freq=agg_dt)

    # Loop over files and pad each to date_min<-->date_max
    for cf in in_files:
        # load it in
        os.chdir(input_dir)
        c_df = pd.read_csv(cf,sep='\t',parse_dates=True) 
        c_df.set_index('datetime',inplace=True)
        c_df.index = pd.to_datetime(c_df.index)
        # Set -9999 to nan
        c_df.replace(-9999,np.NaN,inplace=True)
        # Aggregate to max time step (if needed)
        c_dt  = c_df.index[1]-c_df.index[0]
        if(c_dt!=agg_dt):
            c_df_ag      = c_df.resample(agg_dt,label='left').mean()   
            c_df_ag['p'] = c_df['p'].resample(agg_dt,label='left').sum()
        else:
            c_df_ag = c_df

        # Check we conserve mean
        if(np.sum(c_df.mean()-c_df_ag.mean())!=0):
            print('Warning, differences in mean after aggregation (check if they are significant)')
            print(cf)
            print(c_df.mean()-c_df_ag.mean())
        # pad
        new_df = c_df_ag.reindex(index=new_index,method='pad',fill_value=np.NaN)
        # save to new file
        os.chdir(output_dir)
        new_df.to_csv(cf, sep='\t',index_label='datetime',index_col=0,date_format='%Y%m%dT%H%M%S',na_rep='-9999')


