#! /bin/bash
# example usage: ./step2_preproc.sh 001 02 [B]
# This script will cleanup files to make them BIDS compatible. 
# The goal here is to fix any errors and/or warnings from the BIDS validator. 
set -e #stop immediately if error occurs
# LOAD IN GLOBAL VARIABLES
source globals.sh   
subj=$1
session=$2
acq=$3  # Optional parameter for acquisition name

# ---- TODO 1: Delete scout images from both files and TSV ----
echo "Deleting scout images..."
# Get list of scout files being deleted (to update TSV later)
scout_files=$(find $bids_dir/sub-$subj/ses-$session/anat -name "*scout*" -printf "%P\n")
# Delete the actual files
find $bids_dir/sub-$subj/ses-$session/anat -name "*scout*" -delete

# Delete duplicate runs if you run multiple -- OPTIONAL FOR YOU
echo "Deleting duplicate files..."
dup_files=$(find $bids_dir/sub-$subj -name "*dup*" -printf "%P\n")
find $bids_dir/sub-$subj -name "*dup*" -delete

# ---- TODO 2: Update TSV file to remove deleted files ----
tsv_file="$bids_dir/sub-$subj/ses-$session/sub-${subj}_ses-${session}_scans.tsv"
if [ -f "$tsv_file" ]; then
    echo "Updating scans.tsv file to remove deleted files..."
    # Create a temporary file
    tmp_file=$(mktemp)
    
    # Copy header line
    head -n 1 "$tsv_file" > "$tmp_file"
    
    # Filter out lines containing "scout" or "dup"
    grep -v -e "scout" -e "dup" "$tsv_file" | tail -n +2 >> "$tmp_file"
    
    # Replace original with filtered file
    mv "$tmp_file" "$tsv_file"
    echo "TSV file updated successfully"
else
    echo "Warning: scans.tsv file not found at $tsv_file"
fi

# if you took AP/PA fieldmaps, here's an example on modifying the output to be bids-compatible
# MAKE SURE YOU MODIFY THE FILENAMES TO MATCH YOUR STUDY'S FILENAMES 
# If you want fmriprep to use your fieldmaps for susceptibility distortion correction, 
# you need to tell fmriprep which fieldmaps to use to correct each functional run. 
# To do this, you add an IntendedFor line to fieldmap json files.

# Determine file paths based on whether acquisition parameter is provided
if [ -z "$acq" ]; then
    # No acquisition parameter provided - use original naming scheme
    AP_json="$bids_dir/sub-$subj/ses-${session}/fmap/sub-${subj}_ses-${session}_dir-AP_epi.json"
    PA_json="$bids_dir/sub-$subj/ses-${session}/fmap/sub-${subj}_ses-${session}_dir-PA_epi.json"
else
    # Acquisition parameter provided - include it in the path
    AP_json="$bids_dir/sub-$subj/ses-${session}/fmap/sub-${subj}_ses-${session}_acq-${acq}_dir-AP_epi.json"
    PA_json="$bids_dir/sub-$subj/ses-${session}/fmap/sub-${subj}_ses-${session}_acq-${acq}_dir-PA_epi.json"
fi

# Check if the files exist
if [ ! -f "$AP_json" ] || [ ! -f "$PA_json" ]; then
    echo "Warning: One or both fieldmap files not found:"
    echo "  AP file: $AP_json"
    echo "  PA file: $PA_json"
    echo "Please check your file paths and naming convention."
    exit 1
fi

# ---- TODO 1: Check if IntendedFor field already exists before adding ----
check_intended_for() {
    local json_file="$1"
    if grep -q "IntendedFor" "$json_file"; then
        echo "ERROR: IntendedFor field already exists in $json_file"
        echo "Please check the file and remove any existing IntendedFor fields before running this script."
        exit 1
    fi
}

# Check both files
check_intended_for "$AP_json"
check_intended_for "$PA_json"

# Create the IntendedFor content
beginning='"IntendedFor": ['
runs="" 
for run in {1..16}; do
    run_formatted=$(printf "%02d" $run)
    runs+="\"ses-${session}/func/sub-${subj}_ses-${session}_task-study_run-${run_formatted}_bold.nii.gz\","
done
runs=${runs%,} # Remove the trailing comma
end="],"
insert="${beginning}${runs}${end}"

# insert IntendedFor field after line 35 (i.e., it becomes the new line 36)
echo "Adding IntendedFor field to fieldmap JSON files..."
sed -i "35 a \\ \\ \\ \\ \\ ${insert}" "$AP_json"
sed -i "35 a \\ \\ \\ \\ \\ ${insert}" "$PA_json"

echo "STEP 2 PREPROC COMPLETE"
# and now finally run run_fmriprep.sh!