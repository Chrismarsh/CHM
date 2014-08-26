#!/bin/bash
cppcheck --enable=all `pwd`/src/ --xml 2>err.xml
cppcheck-htmlreport --file=err.xml --report-dir=`pwd`/cppcheck
cd cppcheck
open index.html
