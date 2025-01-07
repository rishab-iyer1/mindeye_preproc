#! /bin/bash

# example usage: ./step2_preproc.sh 001 02

# This script will cleanup files to make them BIDS compatible. 
# The goal here is to fix any errors and/or warnings from the BIDS validator. 
# TODO make sure intendedfor field doesn't get added more than once to the JSON
# TODO update the .tsv files to delete scouts to match the deleted files

set -e #stop immediately if error occurs

# LOAD IN GLOBAL VARIABLES
source globals.sh   

subj=$1
session=$2

# delete scout images
find $bids_dir/sub-$subj/ses-$session/anat -name "*scout*" -delete

# delete duplicate runs if you run multiple -- OPTIONAL FOR YOU
find $bids_dir/sub-$subj -name "*dup*" -delete

# if you took AP/PA fieldmaps, here's an example on modifying the output to be bids-compatible
# MAKE SURE YOU MODIFY THE FILENAMES TO MATCH YOUR STUDY'S FILENAMES 

# If you want fmriprep to use your fieldmaps for susceptibility distortion correction, 
# you need to tell fmriprep which fieldmaps to use to correct each functional run. 

# To do this, you add an IntendedFor line to fieldmap json files. We provide an example below,
# but keep in mind you need to EDIT THIS FOR YOUR SPECIFIC STUDY (e.g., number of runs, task names, etc.)

beginning='"IntendedFor": ['
runs="" 
for run in {1..16}; do
run_formatted=$(printf "%02d" $run)
runs+="\"ses-${session}/func/sub-${subj}_ses-${session}_task-study_run-${run_formatted}_bold.nii.gz\","
done
runs=${runs%,} # Remove the trailing comma
end="],"

insert="${beginning}${runs}${end}"

echo "35 a \\ \ \ \ \ ${insert}"

# insert IntendedFor field after line 35 (i.e., it becomes the new line 36)
sed -i "35 a \\ \ \ \ \ ${insert}" $bids_dir/sub-$subj/ses-${session}/fmap/sub-${subj}_ses-${session}_dir-AP_epi.json
sed -i "35 a \\ \ \ \ \ ${insert}" $bids_dir/sub-$subj/ses-${session}/fmap/sub-${subj}_ses-${session}_dir-PA_epi.json

echo "STEP 2 PREPROC COMPLETE"

# and now finally run run_fmriprep.sh!