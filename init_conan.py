import os

for dir in os.listdir('conan'):

	cmd = 'conan export conan/%s CHM/stable' % dir
	os.system(cmd)