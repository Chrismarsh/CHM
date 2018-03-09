#!/bin/bash

#via https://graphicdesign.stackexchange.com/a/20937
gifsicle -U output.gif `seq -f "#%g" 0 2 545` -O2 -o output_small.gif
