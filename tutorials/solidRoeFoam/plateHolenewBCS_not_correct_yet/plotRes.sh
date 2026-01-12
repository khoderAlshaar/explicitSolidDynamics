#!/bin/bash

# pyFoamPlotWatcher.py --no-default --frequency=60 <logfile>
# pyFoamPlotWatcher.py --no-default --frequency=60 ${1}

pyFoamPlotWatcher.py --write-files --single-data-files-only --hardcopy log.solidRoeFoam
