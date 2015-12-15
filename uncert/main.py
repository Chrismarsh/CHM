import json
import subprocess
import random
#Proof of concept for how to do uncertainty analysis using the command line interface as well as the config files


# example 1, command line interface
# use the harder module and vary a and b
nsim = 4
b = [random.uniform(1, 3) for x in range(nsim)]
c = [random.uniform(0.001, 0.1) for x in range(nsim)]

prj_path = "CHM.config"

for bb,cc in zip(b,c):
    cf1 = "-c config.Harder_precip_phase.const.b:" + str(bb)
    cf2 = "-c config.Harder_precip_phase.const.c:" + str(cc)
    cf3 = "-c output.mesh.base_name:uncert_test/marmot_"+str(bb)+"_"+str(cc)

    subprocess.check_call(['bin/Debug/CHM %s %s %s %s' % (prj_path, cf1, cf2, cf3)],
                          shell=True, cwd="/Users/chris/Documents/PhD/code/CHM")




# example 2, config files
# fname = "../Harder_precip_phase.config"
#
# with open(fname) as f:
#     cfg = json.load(f)
#
# cfg['const']['b'] = 3.14
#
# newfname = "../Harder_precip_phase_test.config"
# with open(newfname,'w') as f:
#    json.dump(cfg,f,indent=4)
#
# with open("../CHM.config") as f:
#     cfg = json.load(f)
#
# cfg['config']['Harder_precip_phase'] = newfname
# cfg['output']['mesh']['base_name']="marmot_harder"
#
# with open("../CHM_harder.config",'w') as f:
#     json.dump(cfg,f,indent=4)
#
