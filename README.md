

Modular unstructured mesh based hydrological model.
![](https://github.com/Chrismarsh/CHM/blob/master/mesh.png)

#Dependencies
* Matlab (optional)
* Armadillo
* boost
* GNU GSL
* Intel TBB
* CGAL
* Paraview (if building the filter) otherwise VTK
* C++11 compliant compiler

#To build:
    cmake .
    make

#To test:
    cmake .
    make check
    make test


#Trouble shooting
##Matlab
### OSX 
* Create a symbolic link from /usr/bin to the matlab install
* ```sudo ln -s /Applications/MATLAB_R2013a.app/bin/matlab /usr/bin/matlab```

###OSX:
## Triangle
Triangle needs to be edited to remove the fpu_control.h header
    $ cat > /usr/include/fpu_control.h

    #define _FPU_SETCW(cw) // nothing
    #define _FPU_GETCW(cw) // nothing
http://stackoverflow.com/a/4766863/410074

###Linux:
Usage of the matlab engine requires installing csh
##Intel compiler
    /usr/lib/armadillo_bits/config.hpp
comment out l. 173

##VTK
Older versions of VTK may have to patch here
http://review.source.kitware.com/#/c/11956/5/Common/Core/vtkMath.h
when building with C++11 

##Google test
Google test can be patched following

http://stackoverflow.com/questions/4655439/print-exception-what-in-google-test

to print the boost::exception diagnostic information

    diff -r /Users/chris/Documents/PhD/code/CHM/tests/gtest/include/gtest/internal/gtest-internal copy.h /Users/chris/Documents/PhD/code/CHM/tests/gtest/include/gtest/internal/gtest-internal.h
    65,66d64
    < #include <boost/exception/all.hpp>
    < 
    1080,1081c1078
    <     catch (boost::exception &e) { \
    <       std::cout << boost::diagnostic_information(e) << std::endl;  \
    ---
    >     catch (...) { \
