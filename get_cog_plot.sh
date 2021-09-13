#!/bin/bash

##for icra cog plot
## yaxis -0.38 ~ 0.025
source /home/utaro/catkin_ws/jaxon_tutorial/devel/setup.bash
cd /home/utaro/cog_correct

cd /home/utaro/cog_correct/fail_plus_JAXON_RED_2021_09-11_16-42_56
datalogger_plotter_with_pyqtgraph.py -f fail_JAXON_RED_2021_09-11_16-42_56 --start 7750 --length 2750 &

cd /home/utaro/cog_correct/fail_nothing_JAXON_RED_2021_09-11_16-52_10
datalogger_plotter_with_pyqtgraph.py -f fail_JAXON_RED_2021_09-11_16-52_10 --start 2500 --length 2750 


