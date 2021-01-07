#!/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Path to data:
WHERE_TO=./plane_layer/test1 
#...........................
#>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Directory of Coral:
#...........................
CORAL_ROOT=/home/ben/dev/coral



#>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Modify stuff below at your own risk :-)
#...........................
if [ -d "$WHERE_TO" ]; then

#Check whether the directory exists
#in which case, cancel the creation

echo " ..."
echo " >>> !!! This directory already exists !!! <<< "
echo " ..."
echo " Aborting the creation of directory/subdirectories"
echo " Delete the directory manually before you can create"
echo " it with this script."                                      


else

# If the directory does not exist, 
# proceed with creation

echo " ..."
echo " Creating the directory:"                                       
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
mkdir ${WHERE_TO}/Snapshots 
        
ln -s ${CORAL_ROOT}/build/periodic_box.exe ${WHERE_TO}/.
cp    ${CORAL_ROOT}/etc/python_scripts/periodic_box/*.py  ${WHERE_TO}/.
cp    ${CORAL_ROOT}/etc/input_files_examples/coral.eqs.salt_fingers  ${WHERE_TO}/coral.equations
cp    ${CORAL_ROOT}/etc/input_files_examples/P3_parameters.in  ${WHERE_TO}/.
cp    ${CORAL_ROOT}/etc/input_files_examples/P3_saving_list.in  ${WHERE_TO}/saving_list.in
 
cd ${WHERE_TO}

fi
