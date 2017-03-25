How to use in parallel mode:
1) Edit config file as per normal settings
2) Config file must hvae the variable config_dir set (this is where individual config files are stored)
3) Run `create_parallel_configs.py parallel_test.py` where parallel_test.py is your config file. This will create one
config file for each variable in config_dir. Also it will create parallel_commands.txt in config_dir, which is the
list of executables that will be past to GNU parallel.
4) Run `run_vtu2geo_in_parallel.sh /config_dir/parallel_commands.txt 10` where 10 is the number of threads to use.

