
prefix=/home/cbm038/libs/

#CGAL
wget https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.7/CGAL-4.7.tar.gz
tar -xf CGAL-4.7.tar.gz
rm CGAL-4.7.tar.gz
cd CGAL-4.7
cmake .
make
cd ..

#boost
http://skylineservers.dl.sourceforge.net/project/boost/boost/1.59.0/boost_1_59_0.tar.gz
tar -xf boost_1_59_0.tar.gz
rm boost_1_59_0.tar.gz
cd boost_1_59_0
./bootstrap
./b2 install --prefix=$prefix/boost159 --layout=versioned
cd ..

#tbb
# wget https://www.threadingbuildingblocks.org/sites/default/files/software_releases/source/tbb44_20151115oss_src.tgz
# tar -xf tbb44_20151115oss_src.tgz
# rm tbb44_20151115oss_src.tgz
# cd tbb44_20151115oss
# make
# cd ..
# mkdir $prefix/tbb44
# mkdir $prefix/tbb44/lib
# cp -r tbb44_20151115oss/include $prefix/tbb44/
# cp -r tbb44_20151115oss/build/`ls | grep *_release`/* $prefix/tbb44/lib

#meteoio
# wget https://models.slf.ch/p/meteoio/downloads/get/MeteoIO-2.5.0-src.tar.gz
# tar -xf MeteoIO-2.5.0-src.tar.gz
# rm MeteoIO-2.5.0-src.tar.gz
# cd MeteoIO-2.5.0-src
# cmake .
# make
# cd ..


