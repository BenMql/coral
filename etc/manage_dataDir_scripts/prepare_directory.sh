#!/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Path to data:
WHERE_TO=./test1 
#...........................
#>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Directory of Coral:
#...........................
CORAL_ROOT=/home/ben/dev/coral

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
mkdir ${WHERE_TO}/CheckPoints
mkdir ${WHERE_TO}/Geometry
mkdir ${WHERE_TO}/Matrices
mkdir ${WHERE_TO}/Profiles
mkdir ${WHERE_TO}/Restart
mkdir ${WHERE_TO}/Slices
mkdir ${WHERE_TO}/Timeseries
mkdir ${WHERE_TO}/Volumes
        
cp ${CORAL_ROOT}/build/coral_LP.exe ${WHERE_TO}/.
cp ${CORAL_ROOT}/etc/python_scripts/*.py ${WHERE_TO}/.
cp ${CORAL_ROOT}/etc/coral_equations_examples/coral.equations.RBC_primitiveVars_Vanilla ${WHERE_TO}/coral.equations
cp ${CORAL_ROOT}/etc/input_files_examples/* ${WHERE_TO}/.
 
cd ${WHERE_TO}

fi
