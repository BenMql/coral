#!/bin/bash

for TARGETFILE in "${LIST_OF_TARGET_FILES[@]}"; do     echo $TARGETFILE && python dealiase_volumes.py $TARGETFILE && rm -f $TARGETFILE/Volumes/*full.dat && rm -f $TARGETFILE/Restart/* && tar -cvf $TARGETFILE.tar $TARGETFILE ; done


