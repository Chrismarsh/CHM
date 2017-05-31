#!/bin/bash
cppcheck --language=c++ --enable=all `pwd`/src/ `pwd`/mesher/ --xml 2>err.xml
./cppcheck-htmlreport.py --file=err.xml --report-dir=`pwd`/cppcheck
cd cppcheck
open index.html
