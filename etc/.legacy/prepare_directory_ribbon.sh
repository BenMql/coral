#!/bin/bash

#WHERE_TO=/home/ben/data/plane_layer/E6R20P07_A_04
#WHERE_TO=/home/ben/data/plane_layer/E6R50P01/E6R50P01_A_05
#WHERE_TO=/home/ben/data/plane_layer/E7R20P07/E7R20P07_B_11
#WHERE_TO=/home/ben/data/plane_layer/E5R20P01_f45/E5R20P01_f45_B_02
#WHERE_TO=/home/ben/data/plane_layer/E7R20P01/E7R20P01_B_09
WHERE_TO=/data/ribbon/R12P1_L2_ell0p1_A02
CORAL_ROOT=/home/bmiquel/dev/coral
echo ${WHERE_TO}

mkdir ${WHERE_TO}
mkdir ${WHERE_TO}/Volumes
mkdir ${WHERE_TO}/XY_Slices
mkdir ${WHERE_TO}/XZ_Slices
mkdir ${WHERE_TO}/YZ_Slices
mkdir ${WHERE_TO}/X_Profiles
mkdir ${WHERE_TO}/Y_Profiles
mkdir ${WHERE_TO}/Z_Profiles
mkdir ${WHERE_TO}/Point_Probes
mkdir ${WHERE_TO}/CheckPoints
mkdir ${WHERE_TO}/Restart
mkdir ${WHERE_TO}/Matrices
mkdir ${WHERE_TO}/Timeseries
        
ln -s ${CORAL_ROOT}/build/ribbon.x ${WHERE_TO}/.
ln -s ${CORAL_ROOT}/etc/python_scripts/*.py  ${WHERE_TO}/.
cp    ${CORAL_ROOT}/etc/input_files_examples/CR_parameters.in  ${WHERE_TO}/.
 
cd ${WHERE_TO}
