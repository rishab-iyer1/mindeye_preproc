#! /bin/bash

# LOAD GLOBAL VARIABLES AND MODULES ON THE CLUSTER
source globals.sh   
# pip install pydeface  # make sure youve pip installed pydeface

sid=$1
session=$2

# # deface T1
T1=`find $bids_dir/sub-$sid/ses-$session/anat -name "*T1w.nii.gz"`
pydeface $T1

# move defaced T1 to extra directory
T1_defaced=`find $bids_dir/sub-$sid/ses-$session/anat -name "*T1w_defaced.nii.gz"`
mv $T1_defaced $defaced_dir/