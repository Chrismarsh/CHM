import os

for dir in os.listdir('conan'):

	cmd = 'conan export conan/%s CHM/dev' % dir
	os.system(cmd)